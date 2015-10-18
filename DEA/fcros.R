#!/usr/bin/Rscript
#
# This script takes an SQLite database obtained from the PIQMIe service and
# performs differential analysis on normalized SILAC protein ratios using
# the FCROS method (Dembele and Kastner, 2013).
#
#
# Author: Arnold Kuzniar
#
# Version: 1.0
#

library(RSQLite, quietly = TRUE)
library(fcros, quietly = TRUE)

fcrosVersion<-packageVersion('fcros')

# check the fcros version
if (fcrosVersion != '1.2') {
   cat('fcros version 1.2 required.\n')
   q(save = 'no', status = 1)
}

scriptname<-basename(sub('.*=', '', commandArgs()[4])) # script name
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
dbfile<-args[1] # SQLite database file
dts<-args[2]  # stores the selected data set (VH10, U2OS or IB10)
fc_cutoff<-as.numeric(args[3]) # robust average or median fold-change (FC2) cutoff (>=1)
p_cutoff<-as.numeric(args[4])  # P-value cutoff (0..1) N.B.: P-values obtaned by this method do not need to be adjusted!
outfile<-file.path(dirname(dbfile), paste0(dts, '_fcros.tab')) # fullpath of output file
nohits<-'Nothing to write onto outfile.\n'

# check user input
if (length(args) != 4 ||
    !dts %in% c('VH10', 'U2OS', 'IB10') ||
    fc_cutoff < 1 ||
    p_cutoff > 1 ||
    p_cutoff < 0) {
    printUsage()
}

# dynamically construct SQL query
# Note: select consistent treat1-4 ratios but do not filter upon ctrl1-2 ratios yet, instead
#       do this post-hoc so that sufficient number of features is available for DE analysis
query<-paste0("SELECT
    A.grp_id grp_id,
    IFNULL(GROUP_CONCAT(DISTINCT gene), '-') genes,
   ", dts, "_L0_M0_H1_norm_ratio_HL AS ratio_H1L0, -- norm. ratio ON/OFF  (treat1)
   ", dts, "_L1_M1_H0_norm_ratio_LH AS ratio_L1H0, -- norm. ratio ON/OFF  (treat2)
   ", dts, "_L0_M0_H1_norm_ratio_HM AS ratio_H1M0, -- norm. ratio ON/OFF  (treat3)
   ", dts, "_L1_M1_H0_norm_ratio_MH AS ratio_M1H0, -- norm. ratio ON/OFF  (treat4)
   ", dts, "_L0_M0_H1_norm_ratio_LM AS ratio_L0M0, -- norm. ratio OFF/OFF (ctrl1)
   ", dts, "_L1_M1_H0_norm_ratio_LM AS ratio_L1M1  -- norm. ratio ON/ON   (ctrl2)
FROM
   VVV_PGROUP_QUANT A, PROT2GRP B, V_PROTEIN C
WHERE
   A.grp_id = B.grp_id
   AND B.prot_acc = C.acc
   AND ((", dts, "_L1_M1_H0_norm_ratio_HL > 1
   AND   ", dts, "_L1_M1_H0_norm_ratio_HM > 1
   AND   ", dts, "_L0_M0_H1_norm_ratio_LH > 1
   AND   ", dts, "_L0_M0_H1_norm_ratio_MH > 1)
   OR   (", dts, "_L1_M1_H0_norm_ratio_LH > 1
   AND   ", dts, "_L1_M1_H0_norm_ratio_MH > 1
   AND   ", dts, "_L0_M0_H1_norm_ratio_HL > 1
   AND   ", dts, "_L0_M0_H1_norm_ratio_HM > 1))
GROUP BY A.grp_id")

# fetch data from db
drv<-dbDriver('SQLite')
conn<-tryCatch(dbConnect(drv, dbname = dbfile),
               error = function(e) {
                  cat(paste0("Error: Database file '", dbfile, "' not found.\n"))
                  q(save = 'no', status = 1)
              })

res<-tryCatch(dbSendQuery(conn, query), error = function(e) {
   cat('Error: Selected data set not found in dbfile.\n')
   q(save = 'no', status = 1)})

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

# columns with log2 ratios corresponding to control conditions
ctrl<-c('log2ratio_L0M0', 'log2ratio_L1M1')

# select columns with log2 ratios corresponding to treatment conditions
treat<-log2cols[!(log2cols %in% ctrl)]

# preliminary raking of proteins based on log2 ratios (treatment conditions only)
FCROS.out<-fcrosMod(tbl[treat], log2.opt = 0, trim.opt = 0.25) 

# change FC2 domain to >=1 so that the magnitudes of up-/down-reg. proteins are on the same scale
FCROS.out$FC2<-ifelse(FCROS.out$FC2 >= 1, FCROS.out$FC2, 1 / FCROS.out$FC2)

# join two data frames and filter post-hoc: log2 ratios treat > ctrl, fold-change, P-values
FCROS.fout<-subset(cbind(tbl, do.call(cbind, FCROS.out[c('FC2', 'ri', 'p.value', 'f.value')])),
   FC2 > fc_cutoff &
   p.value < p_cutoff &
   abs(log2ratio_L1M1) < abs(log2ratio_H1L0) &
   abs(log2ratio_L1M1) < abs(log2ratio_L1H0) &
   abs(log2ratio_L1M1) < abs(log2ratio_H1M0) &
   abs(log2ratio_L1M1) < abs(log2ratio_M1H0) &
   abs(log2ratio_L0M0) < abs(log2ratio_H1L0) &
   abs(log2ratio_L0M0) < abs(log2ratio_L1H0) &
   abs(log2ratio_L0M0) < abs(log2ratio_H1M0) &
   abs(log2ratio_L0M0) < abs(log2ratio_M1H0))

# select columns with numberical values
numcols<-names(FCROS.fout)[!(names(FCROS.fout) %in% skipcols)]

# sort proteins by P-values and round the values in numerical columns 
FCROS.fout<-FCROS.fout[order(FCROS.fout$p.value),]
FCROS.fout[numcols]<-round(FCROS.fout[numcols], 4)

# write results (including header) into tab-delim text file
npgrps.reg<-nrow(FCROS.fout) # count diff. reg. protein groups
cat('Ndiff =', npgrps.reg, '\n')

if (npgrps.reg == 0) {
   warning(nohits)
} else {
   write.table(FCROS.fout, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')
}
exitGracefully(conn, res)
