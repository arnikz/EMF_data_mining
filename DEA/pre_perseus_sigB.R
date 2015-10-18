#!/usr/bin/Rscript
#
# This script is used to pre-process a MaxQuant 'proteinGroups.txt' file
# to streamline the process of adding peak intensity-based significance B
# (sigB) estimates to SILAC ratios using the Perseus software.
#
#
# Author: Arnold Kuzniar
#
# Version: 1.0
#

scriptname<-basename(sub('.*=', '', commandArgs()[4])) # script name
printUsage<-function() {
   cat('Usage:', scriptname, '[proteinGroups.txt]\n')
   q(save = 'no', status = 1)
}

args<-commandArgs(TRUE)
if (length(args) != 1) { printUsage() }

infile<-args[1] # 'proteinGroups.txt' file
outfile<-paste0(infile, '.new')
colnm.prefix<-'Log2 '

tbl<-read.table(infile, header = TRUE, sep = '\t', check.names = FALSE)        # read file into data frame
exp.names<-gsub('Experiment ', '', names(tbl[grep('Experiment', names(tbl))])) # get experiment names
cat(length(exp.names), 'experiments found:', exp.names, '\n')
ratio.cols<-names(tbl[grep('Ratio [HML]/[HML] normalized ', names(tbl))])
log2ratio.cols<-paste0(colnm.prefix, ratio.cols)
intensity.cols<-names(tbl[grep('Intensity [HML] ', names(tbl))])
intensity.H.cols<-intensity.cols[grep('Intensity H ', intensity.cols)]
intensity.M.cols<-intensity.cols[grep('Intensity M ', intensity.cols)]
intensity.L.cols<-intensity.cols[grep('Intensity L ', intensity.cols)]
intensity.HM.cols<-gsub('Intensity H ', 'Intensity H+M ', intensity.H.cols)
intensity.HL.cols<-gsub('Intensity H ', 'Intensity H+L ', intensity.H.cols)
intensity.ML.cols<-gsub('Intensity M ', 'Intensity M+L ', intensity.M.cols)

# add new columns
tbl[log2ratio.cols]<-log2(tbl[ratio.cols])                            # log2 transform ratios
tbl[intensity.HM.cols]<-tbl[intensity.H.cols] + tbl[intensity.M.cols] # sum H and M intensities
tbl[intensity.HL.cols]<-tbl[intensity.H.cols] + tbl[intensity.L.cols] # sum H and L intensities
tbl[intensity.ML.cols]<-tbl[intensity.M.cols] + tbl[intensity.L.cols] # sum M and L intensities

# filter the table for decoys/contaminants and write it onto file
tbl.filtered<-subset(tbl, `Reverse` != '+' & `Contaminant` != '+')
cat(nrow(tbl), 'and', nrow(tbl.filtered), 'protein groups before and after removing decoys/contaminants, respectively.\n')
write.table(tbl.filtered, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, na = 'NaN', sep = '\t')
cat('Done.\n')
