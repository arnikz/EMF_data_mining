#!/usr/bin/Rscript
#
# This script takes an SQLite database obtained from the PIQMIe service and
# performs differential analysis on normalized SILAC protein ratios using
# the RankProd method (Breitling et al., 2004).
#
#
# Author: Arnold Kuzniar
#
# Version: 1.0
#

library(RSQLite, quietly = TRUE)
library(RankProd, quietly = TRUE)

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
seed<-123       # seed for reproducibility
method<-'pfp'   # adjust P-values for multiple hypotheses using the FDR (or Benjamini-Hochberg) method 
dbfile<-args[1] # SQLite database file
dts<-args[2]    # stores the selected data set (VH10, U2OS or IB10)
fc_cutoff<-as.numeric(args[3]) # mean fold-change (FC) cutoff (>=1)
p_cutoff<-as.numeric(args[4])  # adjusted P-value cutoff (0..1)
n_perm<-1000    # number of permutations (default 1000)
outfile<-file.path(dirname(dbfile), paste0(dts, '_rankprod_fc_', fc_cutoff, '_p_', p_cutoff, '.tab')) # fullpath of output file
nohits<-'Nothing to write onto outfile.\n'

# check user input
if (length(args) != 4 ||
    !dts %in% c('VH10', 'U2OS', 'IB10') ||
    fc_cutoff < 1 ||
    p_cutoff > 1 ||
    p_cutoff < 0) {
    printUsage()
}

if (!file.exists(dbfile)) {
   cat(paste0("Error: Database file '", dbfile, "' not found.\n"))
   q(save='no', status=1)
}

# show info
cat('dbfile =', dbfile,
   '\noutfile =', outfile,
   '\ndataset =', dts,
   '\nFC cutoff =', fc_cutoff,
   '\nP-value cutoff =', p_cutoff,
   '\nNperm =', n_perm, '\n')

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
   cat('Error: Selected data set not found.\n')
   q(save = 'no', status = 1) })
tbl<-fetch(res, n = -1)
npgrps.all<-nrow(tbl) # count protein groups
cat('N =', npgrps.all, '\n')

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

# construct class labels of the treated samples
class<-rep(1, length(treat))

# preliminary raking of proteins based on log2 ratios (treatment conditions only)
RP.out<-topGene(RP(tbl[treat], class, logged = TRUE, gene.names = tbl$grp_id, na.rm = FALSE, num.perm = n_perm,
   rand = seed), method = method, num.gene = npgrps.all) # Note: use unique grpIDs instead of gene names

# merge the two tables with up- and down-regulated proteins
RP.df<-as.data.frame(rbind(RP.out$Table1, RP.out$Table2))

# rename fold-change column to 'FC': contains now both class1/class2 and class2/class1 values
names(RP.df)[names(RP.df) == 'FC:(class1/class2)']<-'FC'

# change FC domain to >=1 so that the magnitudes of up-/down-reg. proteins are on the same scale
RP.df$FC<-ifelse(RP.df$FC >= 1, RP.df$FC, 1 / RP.df$FC)

# post-hoc filter: log2 ratios treat > ctrl
df<-subset(tbl,
   abs(log2ratio_L1M1) < abs(log2ratio_H1L0) &
   abs(log2ratio_L1M1) < abs(log2ratio_L1H0) &
   abs(log2ratio_L1M1) < abs(log2ratio_H1M0) &
   abs(log2ratio_L1M1) < abs(log2ratio_M1H0) &
   abs(log2ratio_L0M0) < abs(log2ratio_H1L0) &
   abs(log2ratio_L0M0) < abs(log2ratio_L1H0) &
   abs(log2ratio_L0M0) < abs(log2ratio_H1M0) &
   abs(log2ratio_L0M0) < abs(log2ratio_M1H0))

# add new column 'gene.index', needed to merge tables
df$gene.index<-rownames(df)

# filter by P-value and FC; merge tables on common 'gene.index' column
RP.fout<-merge(df, subset(RP.df, pfp < p_cutoff & FC > fc_cutoff), by='gene.index')

# delete 'gene.index' column
RP.fout$gene.index<-NULL # 'gene.index' column not needed anymore so delete it

# select columns with numberical values
numcols<-names(RP.fout)[!(names(RP.fout) %in% skipcols)]

# sort proteins by adjusted P-values and round the values in numerical columns 
RP.fout<-RP.fout[order(RP.fout$pfp),] # Note: pfp=FDR
RP.fout[numcols]<-round(RP.fout[numcols], 4)

# write results (including header) into tab-delim text file
npgrps.reg<-nrow(RP.fout) # count diff. reg. protein groups
cat('Ndiff =', npgrps.reg, '\n')

if (npgrps.reg == 0) {
   warning(nohits)
} else {
   write.table(RP.fout, file=outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')
}
exitGracefully(conn, res)
