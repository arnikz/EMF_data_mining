#!/usr/bin/env python
#
# This script takes a database (SQLite) obtained from the PIQMIe service and populates
# additional tables/views to facilitate differential protein expression analyses based
# on standardized SILAC ratios.
#
# Note:
#   z_score_{raw|norm}_ratio  - column with canonical Z-score transformed raw/normalized
#                               SILAC protein ratios
#
#   mz_score_{raw|norm}_ratio - column with modified Z-score transformed SILAC protein ratios
#                               suitable for heavy-tailed data (Iglewicz and Hoaglin, 1993)
# 
#
# Author: Arnold Kuzniar
#
# Version: 1.0
#

import os
import sys
import argparse as argp
import math
import numpy as np
import scipy.stats as st
import sqlite3 as sqlt
import collections as cls

ratio_types = { # lookup to link column values to column names
    'RATIO H/L': 'raw_ratio_HL',
    'RATIO H/M': 'raw_ratio_HM',
    'RATIO M/L': 'raw_ratio_ML',
    'RATIO H/L NORMALIZED': 'norm_ratio_HL',
    'RATIO H/M NORMALIZED': 'norm_ratio_HM',
    'RATIO M/L NORMALIZED': 'norm_ratio_ML'
}

score_types = { # lookup to link user input to table column
    'Zr' : 'z_score_raw_ratio',
    'Zn' : 'z_score_norm_ratio',
    'Mr' : 'm_score_raw_ratio',
    'Mn' : 'm_score_norm_ratio'
}

# parse command-line args
parser = argp.ArgumentParser(
   description = 'Differential analysis of SILAC protein ratios based on standardized scores.')

parser.add_argument(
   '-n',
   action = 'store_true',
   dest = 'new_tabs',
   help = 'populate new db tables with (modified) Z-scores')

parser.add_argument(
   '-d',
   dest = 'dts',
   required = True,
   choices = ['VH10', 'U2OS', 'IB10'],
   help = 'select one of the data sets or cell lines')

parser.add_argument(
   '-s',
   required = True,
   choices = score_types.keys(),
   help = 'select one of the score types for filtering: Z*,M* - Z-score or modified Z-score; *r,*n - score based on raw or normalized SILAC protein ratios')

parser.add_argument(
   '-c',
   required = True,
   dest = 'cutoff',
   type = float,
   help = 'absolute score cutoff (e.g. 1.65, 1.96 or 2.58)')

parser.add_argument(
   '-o',
   dest = 'outfile',
   help = 'output file')

parser.add_argument(
   'dbfile',
   help = 'sqlite3 database file')

args = parser.parse_args()

# check user input
dbfile = args.dbfile
outfile = args.outfile
new_tabs = args.new_tabs
dts = args.dts
stype = args.s
cutoff = args.cutoff
score_type = None

if os.path.isfile(dbfile) is False:
   parser.error("dbfile '%s' not found" % dbfile)

if stype and cutoff:
   score_type = score_types[stype]
else:
   parser.error('-s and -c args must be used together')

if outfile is None:
   # set the default output filename
   outfile = os.path.join(os.path.dirname(dbfile), '%s_mzscore_%s_%.2f.tab' % (dts, stype, cutoff))

if cutoff < 0:
   parser.error('the absolute score cutoff must be a positive value')

# print info into STDOUT
print """
dbfile = %s
outfile = %s
dataset = %s
re-score = %s
score type = %s
score cutoff = %.2f
""" % (dbfile, outfile, dts, new_tabs, stype, cutoff)

# sqlite3 user-defined functions (UDFs)
def log(value, base):
    try:
        return math.log(value) / math.log(base)
    except:
        return None

def sqrt(value):
    try:
        return math.sqrt(value)
    except:
        return None

def pvalue(score): # convert Z- or M-score to two-tailed probability (P-value)
    try:
        return 2 * st.norm.cdf(-abs(score))
    except:
        return None

class Stdev: # sample standard deviation (aggregate function)
    def __init__(self):
        self.vec = []

    def step(self, value):
        self.vec.append(value)

    def finalize(self):
        return np.array(self.vec).std(ddof=1)

class Median: # median (aggregate function)
    def __init__(self):
        self.arr = []

    def step(self, value):
        self.arr.append(value)

    def finalize(self):
        return np.median(np.array(self.arr))

class Mad: # median absolute deviation (aggregate function)
    def __init__(self):
        self.arr = []

    def step(self, value):
        self.arr.append(value)

    def finalize(self):
        median = np.median(np.array(self.arr))
        return np.median(np.abs(self.arr - median))

# SQL statements to populate tables/views
sql_create_tables = """
DROP VIEW IF EXISTS V_PGROUP_RATIO;
CREATE VIEW V_PGROUP_RATIO AS
-- simplifies the selection of SILAC ratios/types
SELECT
    A.grp_id,
    exp_name,
    CAST(CASE %s
    END AS TEXT) AS ratio_type,
    CAST(quant_value AS NUMERIC) AS ratio_value
FROM
    PGROUP_QUANT A, V_PGROUP B
WHERE
    A.grp_id = B.grp_id
    AND quant_type IN ('%s')
    AND quant_value;

DROP TABLE IF EXISTS PGROUP_LOG2RATIO_STAT;
CREATE TABLE PGROUP_LOG2RATIO_STAT AS
-- stores descriptive statistics on SILAC protein ratios for each experiment
SELECT
    exp_name,
    ratio_type,
    CAST(COUNT(ratio_value) AS INT) AS n,
    CAST(MIN(LOG(ratio_value, 2)) AS NUMERIC) AS min,
    CAST(MAX(LOG(ratio_value, 2)) AS NUMERIC) AS max,
    CAST(AVG(LOG(ratio_value, 2)) AS NUMERIC) AS mean,
    CAST(MEDIAN(LOG(ratio_value, 2)) AS NUMERIC) AS median,
    CAST(STDEV(LOG(ratio_value, 2)) AS NUMERIC) AS sd,
    CAST(MAD(LOG(ratio_value, 2)) AS NUMERIC) AS mad
FROM
    V_PGROUP_RATIO
GROUP BY
    exp_name, ratio_type;
CREATE INDEX idx_PGROUP_LOG2RATIO_STAT_exp_name_ratio_type ON PGROUP_LOG2RATIO_STAT(exp_name, ratio_type);

DROP VIEW IF EXISTS V_PGROUP_LOG2RATIO_STAT;
CREATE VIEW V_PGROUP_LOG2RATIO_STAT AS
-- shows rounded values of the statistics
SELECT
    exp_name,
    ratio_type,
    n,
    ROUND(min, 4) AS min,
    ROUND(max, 4) AS max,
    ROUND(mean, 4) AS mean,
    ROUND(median, 4) AS median,
    ROUND(sd, 4) AS sd,
    ROUND(mad, 4) AS mad
FROM
    PGROUP_LOG2RATIO_STAT;

DROP TABLE IF EXISTS PGROUP_MZSCORE;
CREATE TABLE PGROUP_MZSCORE AS
-- stores (modified) Z-score transformed SILAC protein raw/norm ratios
SELECT
    grp_id,
    A.exp_name AS exp_name,
    CAST(A.ratio_type AS TEXT) AS ratio_type,
    CAST((LOG(ratio_value, 2) - mean) / sd AS NUMERIC) AS z_score,
    CAST(0.6745 * (LOG(ratio_value, 2) - median) / mad AS NUMERIC) AS m_score
FROM
    V_PGROUP_RATIO A, PGROUP_LOG2RATIO_STAT B
WHERE
    A.exp_name = B.exp_name
    AND A.ratio_type = B.ratio_type;
CREATE INDEX idx_PGROUP_MZSCORE_grp_id ON PGROUP_MZSCORE(grp_id);
CREATE INDEX idx_PGROUP_MZSCORE_exp_name_ratio_type ON PGROUP_MZSCORE(exp_name, ratio_type);
""" % (' '.join([ "\n\tWHEN quant_type='%s' THEN '%s'" % (k, v) for (k, v) in ratio_types.iteritems() ]),
       "','".join(ratio_types.keys()))

# dynamically construct SQL query to select diff. reg. protein groups
sql_sel_pgrps = """
SELECT
    A.grp_id grp_id,
    IFNULL(GROUP_CONCAT(DISTINCT gene), '-') genes,
    {dts}_L0_M0_H1_{score_type}_HL '{stype}_H1L0', -- Z or M-score ON/OFF  (treat1)
    {dts}_L1_M1_H0_{score_type}_LH '{stype}_L1H0', -- Z or M-score ON/OFF  (treat2)
    {dts}_L0_M0_H1_{score_type}_HM '{stype}_H1M0', -- Z or M-score ON/OFF  (treat3)
    {dts}_L1_M1_H0_{score_type}_MH '{stype}_M1H0', -- Z or M-score ON/OFF  (treat4)
    {dts}_L0_M0_H1_{score_type}_LM '{stype}_L0M0', -- Z or M-score OFF/OFF (ctrl1)
    {dts}_L1_M1_H0_{score_type}_LM '{stype}_L1M1', -- Z or M-score ON/ON   (ctrl2)
    PVALUE({dts}_L0_M0_H1_{score_type}_HL) 'pval_H1L0', -- P-value ON/OFF  (treat1)
    PVALUE({dts}_L1_M1_H0_{score_type}_LH) 'pval_L1H0', -- P-value ON/OFF  (treat2)
    PVALUE({dts}_L0_M0_H1_{score_type}_HM) 'pval_H1M0', -- P-value ON/OFF  (treat3)
    PVALUE({dts}_L1_M1_H0_{score_type}_MH) 'pval_M1H0', -- P-value ON/OFF  (treat4)
    PVALUE({dts}_L0_M0_H1_{score_type}_LM) 'pval_L0M0', -- P-value OFF/OFF (ctrl1)
    PVALUE({dts}_L1_M1_H0_{score_type}_LM) 'pval_L1M1'  -- P-value ON/ON   (ctrl2)
FROM
    V_PGROUP_MZSCORE A, PROT2GRP B, V_PROTEIN C
WHERE
    A.grp_id = B.grp_id
    AND B.prot_acc = C.acc
    AND (({dts}_L0_M0_H1_{score_type}_HL > {cutoff}
    AND {dts}_L0_M0_H1_{score_type}_HM > {cutoff}
    AND {dts}_L1_M1_H0_{score_type}_LH > {cutoff}
    AND {dts}_L1_M1_H0_{score_type}_MH > {cutoff})
    OR ({dts}_L0_M0_H1_{score_type}_LH > {cutoff}
    AND {dts}_L0_M0_H1_{score_type}_MH > {cutoff}
    AND {dts}_L1_M1_H0_{score_type}_HL > {cutoff}
    AND {dts}_L1_M1_H0_{score_type}_HM > {cutoff}))
    AND {dts}_L0_M0_H1_{score_type}_ML <= {cutoff}
    AND {dts}_L0_M0_H1_{score_type}_LM <= {cutoff}
    AND {dts}_L1_M1_H0_{score_type}_ML <= {cutoff}
    AND {dts}_L1_M1_H0_{score_type}_LM <= {cutoff}
GROUP BY A.grp_id;
""".format(dts=dts, score_type=score_type, stype=stype, cutoff=cutoff)

# connect to db
with sqlt.connect(args.dbfile) as conn:
   conn.row_factory = sqlt.Row # enable column access by name: row['colnm']
   conn.create_function('log', 2, log)
   conn.create_function('sqrt', 1, sqrt)
   conn.create_function('pvalue', 1, pvalue)
   conn.create_aggregate('stdev', 1, Stdev)
   conn.create_aggregate('median', 1, Median)
   conn.create_aggregate('mad', 1, Mad)
   cur = conn.cursor()

   if new_tabs is True: # populate tables/views only with -n option
      cur.executescript(sql_create_tables)
      cur.execute('SELECT DISTINCT exp_name FROM EXPERIMENT')
      exp_names = [ str(r[0]) for r in cur.fetchall() ]
      cur.execute("SELECT DISTINCT ratio_type FROM PGROUP_LOG2RATIO_STAT")
      ratio_types = [ str(r[0]) for r in cur.fetchall() ]
      n = len(exp_names) * len(ratio_types)
      i = 0
      comma = ','

      # create view for selecting diff. reg. proteins
      sql_create_view = """
DROP VIEW IF EXISTS V_PGROUP_MZSCORE;
CREATE VIEW V_PGROUP_MZSCORE AS
SELECT
    grp_id,
"""
      for e in exp_names:
         for r in ratio_types:
            i += 1
            rr = r[:-2] + r[-2:][::-1] # add inverse ratio: {raw|norm}_ratio_HL changed to *_ratio_HL
            if i == n: comma = ''
            sql_create_view += "\tROUND(CAST(GROUP_CONCAT(CASE WHEN exp_name = '{exp}' AND ratio_type = '{ratio}' THEN z_score ELSE NULL END) AS NUMERIC), 4) AS '{exp}_z_score_{ratio}',\n".format(exp=e, ratio=r)
            sql_create_view += "\tROUND(CAST(GROUP_CONCAT(CASE WHEN exp_name = '{exp}' AND ratio_type = '{ratio}' THEN -1 * z_score ELSE NULL END) AS NUMERIC), 4) AS '{exp}_z_score_{iratio}',\n".format(exp=e, ratio=r, iratio=rr)
            sql_create_view += "\tROUND(CAST(GROUP_CONCAT(CASE WHEN exp_name = '{exp}' AND ratio_type = '{ratio}' THEN m_score ELSE NULL END) AS NUMERIC), 4) AS '{exp}_m_score_{ratio}',\n".format(exp=e, ratio=r)
            sql_create_view += "\tROUND(CAST(GROUP_CONCAT(CASE WHEN exp_name = '{exp}' AND ratio_type = '{ratio}' THEN -1 * m_score ELSE NULL END) AS NUMERIC), 4) AS '{exp}_m_score_{iratio}'{comma}\n".format(exp=e, ratio=r, iratio=rr, comma=comma)
      sql_create_view += "FROM PGROUP_MZSCORE GROUP BY grp_id"
      cur.executescript(sql_create_view)
   
   # write results onto tab-delim file
   if dts is not None:
      sep = '\t' # column separator
      n_pgrps = 0 # count diff. reg. protein groups
      with open(outfile, 'w+') as fout:
         try:
            for drow in [ cls.OrderedDict(xi) for xi in cur.execute(sql_sel_pgrps) ]:
               # first output column names
               if n_pgrps == 0:
                  header = sep.join(drow.keys()) + os.linesep
                  fout.write(header)

               # output remaining rows with column values (grp_id, Z-/M-scores and P-values)
               row = drow.values()
               grp_id = str(drow['grp_id'])
               genes = str(drow['genes'])
               scores = [ str(round(float(x), 4)) for x in row[2:] ]
               srow = grp_id + sep + genes + sep + sep.join(scores) + os.linesep
               fout.write(srow)
               n_pgrps += 1

         except sqlt.OperationalError as e:
            sys.stderr.write('Error: Selected data set not found: %s\n' % e)
            sys.exit(1)

      # remove empty outfile
      if os.path.getsize(outfile) == 0:
         print 'Nothing to write onto outfile.'
         os.remove(outfile)
      else:
         print 'Ndiff =', n_pgrps
