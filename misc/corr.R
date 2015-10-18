#!/usr/bin/Rscript
#
# This batch script reads the EMF database files (*.sqlite) to estimate
# the Pearson's correlation coefficient (r) between reciprocal SILAC ratios,
# makes plots and writes the results in the tab-delimited (*.txt) file.
#
# Column descriptions:
#    [1] emf = EM field type (ELF, UMTS or WiFi)
#    [2] cell_line = cell line (human VH10 or U2OS or mouse IB10)
#    [3] ratio1 = SILAC ratio1 (HL, HM, LM)
#    [4] ratio2 = SILAC ratio2 (LH, MH, ML)
#    [5] N = number of protein groups with both SILAC ratios
#    [6] r = estimated Pearson's correlation coefficient
#    [7] r2 = coefficient of determination
#    [8] CI.L = left (lower) limit of X% confidence interval
#    [9] CI.R = right (upper) limit of X% confidence interval
#    [10] t = Student's t statistic
#    [11] df = degrees of freedom
#    [12] p_value = estimated probability
#
#
# Author: Arnold Kuzniar
#
# Version: 1.0
#

scriptname<-basename(sub('.*=', '', commandArgs()[4])) # script name
printUsage<-function() {
   cat('Usage:', scriptname, '[data dir]\n')
   q(save = 'no', status = 1)
}

args<-commandArgs(TRUE)
if (length(args) != 1){ printUsage() }

library(RSQLite, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(scales, quietly = TRUE)

data_dir<-args[1]
method<-'pearson'
conf_level<-0.95
txt_outfile<-'corr_ratios.txt'
plot_outfile<-'corr_ratios.pdf'
drv<-dbDriver('SQLite')
conn<-NULL
tbl_out<-NULL
emf<-c('ELF', 'UMTS', 'WiFi') # EM field types
org2cell<-list(mouse = c('IB10'), human = c('VH10', 'U2OS')) # lookup for cell lines
ratios<-c('HL', 'HM') # SILAC ratios

# reverse string
revStr<-function (str) {
   return(paste(rev(strsplit(str, split = '')[[1]]), collapse = ''))
}

for (e in emf) {
   for (o in names(org2cell)) {
      base_name<-paste0(toupper(e), '_', o)
      dbfile<-file.path(data_dir, base_name, paste0(base_name, '.sqlite'))
      conn<-dbConnect(drv, dbname = dbfile)
      for (c in unlist(org2cell[o])) {
         for (r1 in ratios) {
            r2<-revStr(r1) # label-swap SILAC ratio
            cat('dbfile: ', dbfile, '\n')
            cat('EMF: ', e, '\n')
            cat('Cell line: ', c, '\n')
            cat('Correlate SILAC ratios: ', paste(r1, 'vs.', r2, sep = ' '), '\n')

            sql<-paste0('
               SELECT
                  CAST(', c, '_L0_M0_H1_norm_ratio_', r1, ' AS REAL) ratio1,
                  CAST(', c, '_L1_M1_H0_norm_ratio_', r2, ' AS REAL) ratio2
               FROM
                  VVV_PGROUP_QUANT
               WHERE ',
                  c, '_L0_M0_H1_norm_ratio_', r1, ' AND ',
                  c, '_L1_M1_H0_norm_ratio_', r2)
            res<-dbSendQuery(conn, sql)
            tbl<-fetch(res, n = -1)
            n<-nrow(tbl)
            cat('Number of quantitated proteins: ', n, '\n\n')

            tbl$ratio1<-log2(tbl$ratio1)
            tbl$ratio2<-log2(tbl$ratio2)
            corr<-cor.test(tbl$ratio1,
                          tbl$ratio2,
                          alternative = 'two.sided',
                          method = method,
                          conf_level = conf_level)
            print(corr)
            CI<-corr$conf.int # confidence interval
            row<-data.frame(emf = e,
               cell_line = c,
               ratio1 = r1,
               ratio2 = r2,
               N = nrow(tbl),
               r = corr$estimate,
               r2 = corr$estimate ^ 2, # coefficient of determination (shared variance)
               CI.L = CI[1],
               CI.R = CI[2],
               t = corr$statistic,
               df = corr$parameter,
               p_value = corr$p.value)
            tbl_out<-rbind(tbl_out, row)
            dbClearResult(res)
         }
      }
   }
}

# write the results table onto a tab-delim file
unlink(txt_outfile)
write.table(tbl_out,
            file = txt_outfile,
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = '\t')

df<-data.frame(dts = paste0(tbl_out$emf, ':', tbl_out$cell_line, ':', tbl_out$ratio1),
               ratio_type = tbl_out$ratio1,
               r = tbl_out$r,
               CI.L = tbl_out$CI.L,
               CI.R = tbl_out$CI.R)

# prepare a dot plot with CI error bars
pdf(plot_outfile)
ggplot(data = df, aes(x = dts, y = r, color = ratio_type)) + geom_point() + geom_errorbar(aes(ymin = CI.L, ymax = CI.R)) + coord_flip() + labs(title = '', x = 'Data set', y = "Pearson's correlation r") + scale_y_continuous(breaks = pretty_breaks(n = 10)) + theme_bw() + theme(legend.position = 'none')

dev.off()
on.exit(dbDisconnect(conn))
