#!/usr/bin/Rscript
#
# This script takes an SQLite database obtained from the PIQMIe service and
# performs differential analysis on normalized SILAC protein ratios using linear
# modeling with empirical Bayes approach as implemented in the 'limma' R package
# (Smyth, 2005; McCarthy and Smyth, 2009).
#
#
# Author: Arnold Kuzniar
#
# Version: 1.0
#

library(RSQLite, quietly = TRUE)
library(limma, quietly = TRUE)
library(statmod, quietly = TRUE)

scriptname<-basename(sub(".*=", "", commandArgs()[4])) # script name
printUsage<-function() {
   cat('Usage:', scriptname, '[dbfile] [dataset={VH10,U2OS,IB10}] [FC cutoff] [P-value cutoff]\n')
   q(save = 'no', status = 1)
}

exitGracefully<-function(conn, res) {
   dbClearResult(res)
   on.exit(dbDisconnect(conn))
   q(save = 'no', status = 0)
}

args<-commandArgs(TRUE)
method<-'fdr'   # adjust P-values for multiple testing using the FDR (or Benjamini-Hochberg) method
dbfile<-args[1] # SQLite database file
dts<-args[2]    # stores the selected data set (VH10, U2OS or IB10)
fc_cutoff<-as.numeric(args[3]) # mean fold-change (FC) cutoff (>=1)
p_cutoff<-as.numeric(args[4])  # adjusted P-value cutoff (0..1)
outfile<-file.path(dirname(dbfile), paste0(dts, '_limma.tab'))
nohits<-'Nothing to write onto outfile.\n'

# check user input
if (length(args) != 4 ||
    !dts %in% c('VH10', 'U2OS', 'IB10') ||
    fc_cutoff < 1 ||
    p_cutoff > 1 ||
    p_cutoff < 0) {
    printUsage()
}

# column aliases->names lookup
colnm<-list(
   'ratio_H1L0' = paste0(dts, '_L0_M0_H1_norm_ratio_HL'),
   'ratio_H0L1' = paste0(dts, '_L1_M1_H0_norm_ratio_HL'),
   'ratio_H1M0' = paste0(dts, '_L0_M0_H1_norm_ratio_HM'),
   'ratio_H0M1' = paste0(dts, '_L1_M1_H0_norm_ratio_HM'),
   'ratio_L0M0' = paste0(dts, '_L0_M0_H1_norm_ratio_LM'),
   'ratio_L1M1' = paste0(dts, '_L1_M1_H0_norm_ratio_LM'))

# construct SQL query: select protein groups having all SILAC ratios
query<-paste0("SELECT
   A.grp_id grp_id,
   IFNULL(GROUP_CONCAT(DISTINCT gene), '-') genes,
   ", colnm[[1]], " AS ", names(colnm)[1], ", -- norm. ratio ON/OFF  (treat1)
   ", colnm[[2]], " AS ", names(colnm)[2], ", -- norm. ratio OFF/ON  (treat2)
   ", colnm[[3]], " AS ", names(colnm)[3], ", -- norm. ratio ON/OFF  (treat3)
   ", colnm[[4]], " AS ", names(colnm)[4], ", -- norm. ratio OFF/ON  (treat4)
   ", colnm[[5]], " AS ", names(colnm)[5], ", -- norm. ratio OFF/OFF (ctrl1)
   ", colnm[[6]], " AS ", names(colnm)[6], "  -- norm. ratio ON/ON   (ctrl2)
FROM
   VVV_PGROUP_QUANT A, PROT2GRP B, V_PROTEIN C
WHERE
   A.grp_id = B.grp_id
   AND B.prot_acc = C.acc
   AND ", colnm[[1]], "
   AND ", colnm[[2]], "
   AND ", colnm[[3]], "
   AND ", colnm[[4]], "
GROUP BY A.grp_id")
#cat('SQL query = \n', query, '\n')

# construct targets list: use direct two-color design
# take into account all treatments and controls (including label-swaps)
# Note: direct correspondance of columns selected by the SQL query and Cy3/Cy5 vectors
targets<-cbind( # "squeeze" triplex SILAC to two-color marray desing
   Cy3=c('OFF', 'ON', 'OFF', 'ON', 'OFF', 'ON'),  # Cy3 -> Light or Medium SILAC labels
   Cy5=c('ON', 'OFF', 'ON', 'OFF', 'OFF', 'ON'))  # Cy5 -> Heavy or Light SILAC label
rownames(targets)<-unlist(colnm, use.names = FALSE) # add experiment names including selected SILAC ratios

# construct design/contrast matrix from targets: the resulting coeficients measure the changes relative to ref
design<-modelMatrix(targets, ref = 'OFF', verbose = FALSE)

# fetch data from db
drv<-dbDriver('SQLite')
conn<-tryCatch(dbConnect(drv, dbname = dbfile),
               error = function(e) {
                  cat(paste0("Error: Database file '", dbfile, "' not found.\n"))
                  q(save = 'no', status = 1)
              })

res<-tryCatch(dbSendQuery(conn, query),
              error = function(e) {
                 cat('Error: Selected data set not found in dbfile.\n')
                 q(save = 'no', status = 1)
              })
tbl<-fetch(res, n = -1)
npgrps.all<-nrow(tbl) # count protein groups

# show info
cat('dbfile =', dbfile,
    '\noutfile =', outfile,
    '\ndataset =', dts,
    '\nFC cutoff =', fc_cutoff,
    '\nP-value cutoff =', p_cutoff,
    '\nN =', npgrps.all, '\n')

# exit gracefully if no data returned
if (npgrps.all == 0) {
   warning(nohits)
   exitGracefully(conn, res)
}

#
# DE analysis
#

# columns excluded from further processing
skipcols<-c('grp_id', 'genes')

# dynamically generate column names with log2ratio_* prefix
rcols<-names(tbl)[!(names(tbl) %in% skipcols)] # select columns with SILAC ratios
log2cols<-paste0('log2', rcols)

# log2-transform SILAC ratios, append new columns to the table
tbl[log2cols]<-log2(tbl[rcols])

# take into account different LC-MS runs: two runs each yielding three SILAC ratios
msruns<-rep(c(1, 2), length(colnm) / 2)

# get indices of columns with log2 ratios
log2cols<-grep('log2', names(tbl))
corfit<-duplicateCorrelation(tbl[log2cols], design, block = msruns)

# fit linear model
fit<-lmFit(tbl[log2cols], design, block = msruns, cor = corfit$consensus)

# compute empirical Bayes moderated-t p-values
#fit.eBayes<-eBayes(fit)
# compute empirical Bayes moderated-t p-values relative to a minimum required fold-change threshold (TREAT)
fit.treat<-treat(fit)

# filter the protein list by FDR and log2FC
#LM.df<-topTable(fit.eBayes,
LM.df<-topTreat(fit.treat,
   genelist = tbl$grp_id,
   number = npgrps.all,
   adjust.method = method,
   p.value = p_cutoff,
   lfc = log2(fc_cutoff),
   confint = TRUE)

npgrps.reg<-nrow(LM.df) # count diff. reg. protein groups
cat('Ndiff =', npgrps.reg, '\n')

# exit upon empty list
if (npgrps.reg == 0) {
   warning(nohits)
   exitGracefully(conn, res)
}

# merge tables
names(LM.df)[names(LM.df) == 'ID']<-'grp_id' # rename column
LM.fout<-merge(tbl, LM.df, by = 'grp_id')

# sort proteins by adjusted P-values and round the values in numerical columns 
numcols<-names(LM.fout)[!(names(LM.fout) %in% skipcols)]
LM.fout<-LM.fout[order(LM.fout$adj.P.Val),]
LM.fout[numcols]<-round(LM.fout[numcols], 4)

# write results (including header) into tab-delim text file
write.table(LM.fout, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')
exitGracefully(conn, res)
