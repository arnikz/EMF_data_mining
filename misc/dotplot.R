#!/usr/bin/Rscript
#
# This batch script reads the EMF database files (*.sqlite) and outputs dotplots of
# normalized SILAC protein ratios for control and treatment conditions.
#
#
# Author: Arnold Kuzniar
#
# Version: 1.0
#

scriptname<-basename(sub(".*=", "", commandArgs()[4])) # script name
printUsage<-function() {
   cat('Usage:', scriptname, '[data dir]\n')
   q(save = 'no', status = 1)
}

args<-commandArgs(TRUE)
if (length(args) != 1){ printUsage() }

library(RSQLite)
library(lattice)

data_dir<-args[1]
outfile<-'EMF_dotplot_ratios.pdf'
drv<-dbDriver('SQLite')
conn<-NULL
df<-NULL
emf<-c('WiFi', 'UMTS', 'ELF') # EM field types
org2cell<-list(mouse = c('IB10'), human = c('U2OS', 'VH10')) # lookup table for cell lines
n_rnd_rows = 1000 # number of randonly selected rows

# fetch quantitative data from the database files into a single table
for (e in emf) {
   for (o in names(org2cell)) {
      base_name<-paste0(e, '_', o)
      dbfile<-file.path(data_dir, base_name, paste0(base_name, '.sqlite'))
      conn<-dbConnect(drv, dbname = dbfile)
      for (c in unlist(org2cell[o])) {
         # construct SQL query
         not_null_ratios<-paste0( # part of the WHERE clause
            c, "_L0_M0_H1_norm_ratio_HL AND ",
            c, "_L0_M0_H1_norm_ratio_HM AND ",
            c, "_L0_M0_H1_norm_ratio_LM AND ",
            c, "_L1_M1_H0_norm_ratio_LH AND ",
            c, "_L1_M1_H0_norm_ratio_MH AND ",
            c, "_L1_M1_H0_norm_ratio_LM")

          sql<-paste0("
             SELECT * FROM (
                SELECT
                   'treated' ratio_class,
                   ", c, "_L0_M0_H1_norm_ratio_HL ratio_value
                FROM
                   VVV_PGROUP_QUANT
                WHERE ", not_null_ratios, "
                ORDER BY RANDOM() LIMIT ", n_rnd_rows, " / 8)
             UNION
             SELECT * FROM (
                SELECT
                   'treated' ratio_class,
                   ", c, "_L0_M0_H1_norm_ratio_LH ratio_value
                FROM
                   VVV_PGROUP_QUANT
                WHERE ", not_null_ratios, "
                ORDER BY RANDOM() LIMIT ", n_rnd_rows, " / 8)
             UNION
             SELECT * FROM (
                SELECT
                   'treated' ratio_class,
                   ", c, "_L1_M1_H0_norm_ratio_LH ratio_value
                FROM
                   VVV_PGROUP_QUANT
                WHERE ", not_null_ratios, "
                ORDER BY RANDOM() LIMIT ", n_rnd_rows, " / 8)
             UNION
             SELECT * FROM (
                SELECT
                   'treated' ratio_class,
                   ", c, "_L1_M1_H0_norm_ratio_HL ratio_value
                FROM
                   VVV_PGROUP_QUANT
                WHERE ", not_null_ratios, "
                ORDER BY RANDOM() LIMIT ", n_rnd_rows, " / 8)
             UNION
             SELECT * FROM (
                SELECT
                   'treated' ratio_class,
                   ", c, "_L0_M0_H1_norm_ratio_HM ratio_value
                FROM
                   VVV_PGROUP_QUANT
                WHERE ", not_null_ratios, "
                ORDER BY RANDOM() LIMIT ", n_rnd_rows, " / 8)
             UNION
             SELECT * FROM (
                SELECT
                   'treated' ratio_class,
                   ", c, "_L0_M0_H1_norm_ratio_MH ratio_value
                FROM
                   VVV_PGROUP_QUANT
                WHERE ", not_null_ratios, "
                ORDER BY RANDOM() LIMIT ", n_rnd_rows, " / 8)
             UNION
             SELECT * FROM (
                SELECT
                   'treated' ratio_class,
                   ", c, "_L1_M1_H0_norm_ratio_MH ratio_value
                FROM
                   VVV_PGROUP_QUANT
                WHERE ", not_null_ratios, "
                ORDER BY RANDOM() LIMIT ", n_rnd_rows, " / 8)
             UNION
             SELECT * FROM (
                SELECT
                   'treated' ratio_class,
                   ", c, "_L1_M1_H0_norm_ratio_HM ratio_value
                FROM
                   VVV_PGROUP_QUANT
                WHERE ", not_null_ratios, "
                ORDER BY RANDOM() LIMIT ", n_rnd_rows, " / 8)
             UNION
             SELECT * FROM (
                SELECT
                   'control' ratio_class,
                   ", c, "_L0_M0_H1_norm_ratio_LM ratio_value
                FROM
                   VVV_PGROUP_QUANT
                WHERE ", not_null_ratios, "
                ORDER BY RANDOM() LIMIT ", n_rnd_rows, " / 4)
             UNION
             SELECT * FROM (
                SELECT
                   'control' ratio_class,
                   ", c, "_L0_M0_H1_norm_ratio_ML ratio_value
                FROM
                   VVV_PGROUP_QUANT
                WHERE ", not_null_ratios, "
                ORDER BY RANDOM() LIMIT ", n_rnd_rows, " / 4)
             UNION
             SELECT * FROM (
                SELECT
                   'control' ratio_class,
                   ", c, "_L1_M1_H0_norm_ratio_LM ratio_value
                FROM
                   VVV_PGROUP_QUANT
                WHERE ", not_null_ratios, "
                ORDER BY RANDOM() LIMIT ", n_rnd_rows, " / 4)
             UNION
             SELECT * FROM (
                SELECT
                   'control' ratio_class,
                   ", c, "_L1_M1_H0_norm_ratio_ML ratio_value
                FROM
                   VVV_PGROUP_QUANT
                WHERE ", not_null_ratios, "
                ORDER BY RANDOM() LIMIT ", n_rnd_rows, " / 4)")

         #cat(sql, "\n")
         cat(e, '\t', c, '\n') # print progress info to STDOUT
         res<-dbSendQuery(conn, sql)
         tbl<-fetch(res, n = -1) # fetch rows
  
         df<-rbind(df, # append rows to data frame
            data.frame(
               data_set = paste0(e, ':', c, ':', tbl$ratio_class), # EMF:cell line:SILAC ratio class
               log2_ratio = log2(tbl$ratio_value)))                # log2-transformed SILAC protein ratio
         dbClearResult(res) # clean-up
      }
   }
}
pdf(outfile)
stripplot(data_set ~ log2_ratio,
   df,
   jitter = TRUE,
   col = 'black',
   alpha = .06,
   pch = 20,
   cex = .5,
   grid = 'v',
   main = 'Distribution of SILAC protein ratios',
   xlab = expression(paste(log[2], '(ratio)')),
   ylab = 'Data set',
   xlim = c(-2, 2), # limit log2 ratios
   scales = list(tick.number = 10))
dev.off()
on.exit(dbDisconnect(conn))
