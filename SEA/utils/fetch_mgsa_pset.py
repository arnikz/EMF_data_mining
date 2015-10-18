#!/usr/bin/env python
#
# This script fetches a population set of protein accessions from an SQLite database and
# writes the result set into an output file used for MGSA analysis (Bauer et al., 2010).
#
# Author: Arnold Kuzniar
#
# Version 1.0
#

import os
import sys
import argparse as argp
import sqlite3 as sqlt

def main():
   parser = argp.ArgumentParser(
      description = 'Fetch a population set for the MGSA analysis.')

   parser.add_argument(
      '-d', '--dataset',
      dest = 'dts',
      help = 'select one of the cell lines as input dataset',
      required = True,
      choices = ['VH10', 'U2OS', 'IB10'])

   parser.add_argument(
      '-o', '--outfile',
      dest = 'outfile',
      help = 'output text file e.g. {EMF}_{species|cell_line}_all_acc.txt',
      required = True)

   parser.add_argument(
      'dbfile',
      help = 'path to a dbfile in sqlite3 format')

   args = parser.parse_args()

   if os.path.isfile(args.dbfile) is False:
      parser.error("dbfile '%s' not found" % args.dbfile)

   dts = args.dts
   outfile = args.outfile
   dbfile = args.dbfile

   # Select all protein groups associated with no missing SILAC ratios (six in total).
   # N.B.: Only a leading protein (accession) per group is selected.
   sql_query = """
SELECT
   DISTINCT(prot_acc) AS prot_acc
FROM
   VVV_PGROUP_QUANT INNER JOIN PROT2GRP USING(grp_id) INNER JOIN PEP2PROT USING(prot_acc)
WHERE
   {exp1}_{HL} AND
   {exp1}_{HM} AND
   {exp1}_{ML} AND
   {exp2}_{HL} AND
   {exp2}_{HM} AND
   {exp2}_{ML} AND
   lead_prot != 0;
   """.format(exp1 = dts + '_L0_M0_H1',
              exp2 = dts + '_L1_M1_H0',
              HL   = 'norm_ratio_HL',
              HM   = 'norm_ratio_HM',
              ML   = 'norm_ratio_ML')

   # connect to db and write result set onto file
   with open(outfile, 'w') as fout:
      with sqlt.connect(dbfile) as conn:
         try:
            cur = conn.cursor()
            for row in cur.execute(sql_query):
               fout.write(str(row[0]) + os.linesep)
            cur.close()

         except sqlt.OperationalError as err:
            sys.stderr.write('Error: Selected data set not found in %s: %s\n' % (dbfile, err))
            os.remove(outfile)
            sys.exit(1)

if __name__ == '__main__':
   main()
