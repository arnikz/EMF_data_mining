#!/usr/bin/env python
#
# This script takes a database file (SQLite) obtained from the PIQMIe service to perform differential
# protein expression analysis using fold-changes (normalized SILAC ratios) and peak-intensity based
# significance B (Cox and Mann, 2008).
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
import sqlite3 as sqlt
import collections as cls

def main():
   # parse command-line args
   parser = argp.ArgumentParser(
      description = 'Differential analysis based on normalized SILAC protein ratios and peak-intensity based significance B.')

   parser.add_argument(
      '-d',
      dest     = 'dts',
      required = True,
      choices  = ['VH10', 'U2OS', 'IB10'],
      help     = 'select one of the cell lines as input dataset')

   parser.add_argument(
      '-f',
      dest    = 'fc_cutoff',
      type    = float,
      default = 1,
      help    = 'min. fold-change required to filter SILAC protein ratios')

   parser.add_argument(
      '-s',
      dest    = 'sigB_cutoff',
      type    = float,
      default = 1,
      help    = 'peak-intensity based significance B cutoff (adjusted P-value or FDR)')

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
   dts = args.dts
   fc_cutoff = args.fc_cutoff
   sigB_cutoff = args.sigB_cutoff

   if os.path.isfile(dbfile) is False:
      parser.error("dbfile '%s' not found" % dbfile)

   if fc_cutoff < 1:
      parser.error('fold-change cutoff cannot be smaller than 1')

   if sigB_cutoff < 0 or sigB_cutoff > 1:
      parser.error('sigB cutoff (P-value) must be between 0 and 1')

   if outfile is None:
      outfile = os.path.join(os.path.dirname(dbfile), '%s_fcsigB.tab' % dts) # fullpath of the outfile

   print """
dbfile = %s
outfile = %s
dataset = %s
FC cutoff = %.4f
P-value cutoff = %.4f
   """ % (dbfile, outfile, dts, fc_cutoff, sigB_cutoff)

   # sqlite3 user-defined function (UDF)
   def log(value, base):
      try:
         return math.log(value) / math.log(base)
      except:
         return None

   # construct SQL query to select diff. reg. protein groups
   # Note: filter by fold-change only if sigB cutoff set to 1 (default)

   sigB_filter = """
      AND {dts}_L0_M0_H1_sig_ratio_HL < {cutoff}
      AND {dts}_L0_M0_H1_sig_ratio_HM < {cutoff}
      AND {dts}_L1_M1_H0_sig_ratio_LH < {cutoff}
      AND {dts}_L1_M1_H0_sig_ratio_MH < {cutoff}
      AND {dts}_L0_M0_H1_sig_ratio_LM >= {cutoff}
      AND {dts}_L1_M1_H0_sig_ratio_LM >= {cutoff}
   """.format(dts = dts, cutoff = sigB_cutoff)

   treat_gt_ctrl_ratio_filter = """
      AND MIN(MAX({dts}_L0_M0_H1_norm_ratio_HL, {dts}_L0_M0_H1_norm_ratio_LH),
              MAX({dts}_L0_M0_H1_norm_ratio_HM, {dts}_L0_M0_H1_norm_ratio_MH),
              MAX({dts}_L1_M1_H0_norm_ratio_LH, {dts}_L1_M1_H0_norm_ratio_HL),
              MAX({dts}_L1_M1_H0_norm_ratio_MH, {dts}_L1_M1_H0_norm_ratio_HM)) >
          MAX({dts}_L0_M0_H1_norm_ratio_LM, {dts}_L0_M0_H1_norm_ratio_ML,
              {dts}_L1_M1_H0_norm_ratio_LM, {dts}_L1_M1_H0_norm_ratio_ML)
   """.format(dts = dts)

   extra_filter = sigB_filter if sigB_cutoff != 1 else treat_gt_ctrl_ratio_filter

   sql_sel_pgrps = """
      SELECT
         A.grp_id grp_id,
         IFNULL(GROUP_CONCAT(DISTINCT gene), '-') genes,
         {dts}_L0_M0_H1_norm_ratio_HL AS ratio_H1L0, -- norm. ratio ON/OFF  (treat1)
         {dts}_L1_M1_H0_norm_ratio_LH AS ratio_L1H0, -- norm. ratio ON/OFF  (treat2)
         {dts}_L0_M0_H1_norm_ratio_HM AS ratio_H1M0, -- norm. ratio ON/OFF  (treat3)
         {dts}_L1_M1_H0_norm_ratio_MH AS ratio_M1H0, -- norm. ratio ON/OFF  (treat4)
         {dts}_L0_M0_H1_norm_ratio_LM AS ratio_L0M0, -- norm. ratio OFF/OFF (ctrl1)
         {dts}_L1_M1_H0_norm_ratio_LM AS ratio_L1M1, -- norm. ratio ON/ON   (ctrl2)
         LOG({dts}_L0_M0_H1_norm_ratio_HL, 2) AS log2ratio_H1L0, -- log2 ratio ON/OFF  (treat1)
         LOG({dts}_L1_M1_H0_norm_ratio_LH, 2) AS log2ratio_L1H0, -- log2 ratio ON/OFF  (treat2)
         LOG({dts}_L0_M0_H1_norm_ratio_HM, 2) AS log2ratio_H1M0, -- log2 ratio ON/OFF  (treat3)
         LOG({dts}_L1_M1_H0_norm_ratio_MH, 2) AS log2ratio_M1H0, -- log2 ratio ON/OFF  (treat4)
         LOG({dts}_L0_M0_H1_norm_ratio_LM, 2) AS log2ratio_L0M0, -- log2 ratio OFF/OFF (ctrl1)
         LOG({dts}_L1_M1_H0_norm_ratio_LM, 2) AS log2ratio_L1M1, -- log2 ratio ON/ON   (ctrl2)
         {dts}_L0_M0_H1_sig_ratio_HL AS pval_H1L0, -- sigB ON/OFF  (treat1)
         {dts}_L1_M1_H0_sig_ratio_LH AS pval_L1H0, -- sigB ON/OFF  (treat2)
         {dts}_L0_M0_H1_sig_ratio_HM AS pval_H1M0, -- sigB ON/OFF  (treat3)
         {dts}_L1_M1_H0_sig_ratio_MH AS pval_M1H0, -- sigB ON/OFF  (treat4)
         {dts}_L0_M0_H1_sig_ratio_LM AS pval_L0M0, -- sigB OFF/OFF (ctrl1)
         {dts}_L1_M1_H0_sig_ratio_LM AS pval_L1M1  -- sigB ON/ON   (ctrl2)
     FROM
        VVV_PGROUP_QUANT A, PROT2GRP B, V_PROTEIN C
     WHERE
        A.grp_id = B.grp_id
        AND B.prot_acc = C.acc
        AND (({dts}_L0_M0_H1_norm_ratio_HL > {cutoff}
        AND {dts}_L0_M0_H1_norm_ratio_HM > {cutoff}
        AND {dts}_L1_M1_H0_norm_ratio_LH > {cutoff}
        AND {dts}_L1_M1_H0_norm_ratio_MH > {cutoff})
        OR ({dts}_L0_M0_H1_norm_ratio_LH > {cutoff}
        AND {dts}_L0_M0_H1_norm_ratio_MH > {cutoff}
        AND {dts}_L1_M1_H0_norm_ratio_HL > {cutoff}
        AND {dts}_L1_M1_H0_norm_ratio_HM > {cutoff})) {filter}
     GROUP BY A.grp_id;
   """.format(dts = dts, cutoff = fc_cutoff, filter = extra_filter)

   # connect to db and write result set into file
   with open(outfile, 'w') as fout:
      with sqlt.connect(dbfile) as conn:
         sep = '\t'                  # column separator
         n_pgrps = 0                 # count differentially regulated proteins (groups)
         conn.row_factory = sqlt.Row # enable column access by name: row['colnm']
         conn.create_function('log', 2, log) # register log() UDF

         try:
            cur = conn.cursor()
            for drow in [ cls.OrderedDict(xi) for xi in cur.execute(sql_sel_pgrps) ]:
               # header line with column names
               if n_pgrps == 0:
                  header = sep.join(drow.keys()) + os.linesep
                  fout.write(header)

               # rows with column values: grp_id, log2 ratios and sigB
               row = drow.values()
               grp_id = str(drow['grp_id'])
               genes = str(drow['genes'])
               scores = sep.join([ 'NA' if x is None else str(round(float(x), 4)) for x in row[2:] ])
               srow = sep.join([grp_id, genes, scores]) + os.linesep
               fout.write(srow)
               n_pgrps += 1

            cur.close()

         except sqlt.OperationalError as err:
            sys.stderr.write('Error: Selected data set not found: %s\n' % err)
            os.remove(outfile)
            sys.exit(1)

   # remove empty outfile
   if os.path.getsize(outfile) == 0:
      print 'Nothing to write onto outfile.'
      os.remove(outfile)
   else:
      print 'Ndiff =', n_pgrps

if __name__ == '__main__':
   main()
