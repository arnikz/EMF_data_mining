#!/usr/bin/Rscript
#
# This script performs model-based protein set analysis (Bauer et al., 2010, 2011) using
# pathway/biological process associations from UniProt-GOA, KEGG and Reactome.
#
# Input:
#    population set = all quantitated proteins detected in an EMF exposure experiment
#    study set = differentially regulated proteins detected upon EMF exposure
# 
# Output: 
#    significant categories written onto STDOUT
#    tab-delimited text file containing the probability estimates for all mapped categories
#    PDF file with several plots on MGSA fits
#
#
# Author: Arnold Kuzniar
#
# Version: 1.0
#

library(mgsa, quietly = TRUE)

scriptname<-basename(sub('.*=', '', commandArgs()[4])) # script name
printUsage<-function() {
   cat('Usage:', scriptname, '[db={goa,kegg,reactome}] [marginal posterior prob. cutoff]\n')
   q(save = 'no', status = 1)
}

checkFile<-function(file) {
   if (!file.exists(file)) {
      cat('Error: File', file, 'not found!\n')
      q(save = 'no', status = 1)
   }
}

args<-commandArgs(TRUE)
db<-args[1] # source database
cutoff<-as.numeric(args[2]) # marginal posterior probability cutoff

# check user input
if (length(args) != 2 ||
    cutoff < 0 ||
    cutoff > 1 ||
    !db %in% c('goa', 'kegg', 'reactome')) {
    printUsage()
}

# read protein-to-KEGG/Reactome pathways association file
readPathways<-function (filename) {
   # read tab-delim text file with three columns: [path_id] [path_name] [uniprot_acc]
   df<-read.table(filename, header = TRUE, sep = '\t', quote = '')
   n<-nrow(df)
   sets<-list()
   df$uniprot_acc<-by(
      df,
      1:n,
      # split comma separated UniProtKB protein accessions
      function (row) {
         unlist(strsplit(as.character(row$uniprot_acc), ',', fixed = TRUE))
      }
   )
   for (i in 1:n) {
      # append accession sets to the list
      sets<-c(sets, setNames(list(df$uniprot_acc[[i]]), df$path_id[i]))
      #print(sets)
   }
   mgsa<-new('MgsaSets', sets = sets)
   # add KEGG pathway annotations
   mgsa@setAnnotations<-data.frame(
      row.names = as.character(df$path_id),
      term = as.character(df$path_name))
   return(mgsa)
}

emf<-c('ELF', 'UMTS', 'WIFI', 'RF', 'EMF') # N.B.: RF = pooled UMTS and WIFI data; EMF = pooled ELF, UMTS and WIFI data
org<-c('human', 'mouse')
S.file.sfx<-'diff_acc.txt'
P.file.sfx<-'all_acc.txt'
aspect = 'P' # select GO sub-ontology:
             # P = biological process (default)
             # F = biological function
             # C = cellular component

pdf.outfile<-paste0('mgsa_', db, '.pdf')
txt.outfile<-paste0('mgsa_', db, '.txt')

# delete outfiles from previous run
unlink(pdf.outfile)
unlink(txt.outfile)

pdf(pdf.outfile)
for (o in org) {
   # read ontological categories from association file and instantiate MgsaSets
   gas.file<-file.path(db, paste0('gene_association.', db, '_', o))
   checkFile(gas.file)
   mgsa<-NULL

   if (db == 'goa') {
      mgsa<-readGAF(gas.file, aspect = aspect)
   } else {
      mgsa<-readPathways(gas.file)
   }

   for (e in emf) {
      # show data set info
      cat(toupper(db), ':', e, o, '\n')

      # read study (S) and population (P) sets from file
      S.file<-file.path(db, paste(e, o, S.file.sfx, sep = '_'))
      P.file<-file.path(db, paste(e, o, P.file.sfx, sep = '_'))

      for (f in c(S.file, P.file)) { checkFile(f) }

      S.set<-c(t(read.table(S.file)))
      P.set<-c(t(read.table(P.file)))

      # fit the MGSA model for the input data set
      fit<-mgsa(S.set, mgsa, population = P.set, restarts = 20)

      # provide a graphical representation of the fit for each EMF data set
      plot(fit)
      title(main = paste0(toupper(db), ': ', e, ' ', o), outer = TRUE)

      # show significant categories in STDOUT given the posterior probability cutoff
      res<-setsResults(fit)
      res<-res[order(res[, c('estimate')], decreasing = TRUE),] # sort results by posterior probability
      print(subset(res, estimate > cutoff))

      # write (append) results onto a tab-delimited file
      write.table(data.frame(emf = e, org = o, term_id = row.names(res), subset(res, select = c('estimate', 'std.error', 'inStudySet', 'inPopulation', 'term'))),
            file = txt.outfile,
            append = TRUE,
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = '\t')
   }
}
warnings()
dev.off()
