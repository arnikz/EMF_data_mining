#!/usr/bin/Rscript
#
# This script takes the original MaxQuant 'proteinGroups.txt' file and its
# post-processed version with peak intensity-based significance B (sigB)
# estimates (obtained by the Perseus software) and outputs a single merged
# file suitable for downstream analysis using the PIQMIe service.
#
#
# Author: Arnold Kuzniar
#
# Version: 1.0
#

scriptname<-basename(sub('.*=', '', commandArgs()[4])) # script name
printUsage<-function() {
   cat('Usage:', scriptname, '[FILE1] [FILE2]\n')
   cat(" FILE1 - original 'proteinGroups.txt' file without sigB estimates\n")
   cat(" FILE2 - post-processed 'proteinGroups.txt' file with sigB estimates\n\n")
   q(save = 'no', status = 1)
}

args<-commandArgs(TRUE)
file1<-args[1]
file2<-args[2]

if (!file.exists(file1) || !file.exists(file2)) { printUsage() }
colnm.prefix<-'Log2 '
dir<-dirname(file1)
outfile<-file.path(dir, paste(basename(file1), basename(file2), sep = '-'))
tbl1<-read.table(file1, header = TRUE, sep = '\t', check.names = FALSE) # read file into data frame
tbl2<-read.table(file2, header = TRUE, sep = '\t', check.names = FALSE)
sigB.cols<-names(tbl2[grep('Significance B', names(tbl2))]) # get sigB column names
new.sigB.cols<-gsub(colnm.prefix, '', sigB.cols)            # new sigB column names
names(tbl2)[names(tbl2) %in% sigB.cols]<-new.sigB.cols      # rename sigB columns
join.col<-'Protein IDs'
write.table(merge(tbl1, tbl2[c(join.col, new.sigB.cols)], by = join.col, all.x = TRUE ), # left-outer join
            file = outfile, row.names = FALSE, col.names = TRUE, na = 'NaN',
            quote = FALSE, sep = '\t')
