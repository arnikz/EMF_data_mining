#!/usr/bin/Rscript
#
# This batch script summarizes the variability of protein quantitations for all EMF data sets.
#
# Input:
#    base directory with database files (*.sqlite)
#
# Output:
#    stat_quant.pdf file with bar charts
#    stat_quant.txt file with tabulated statistical data
#    prot_quant.txt file with filtered set of protein quantitations
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

library(RSQLite, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(lattice, quietly = TRUE)

data_dir<-args[1]
plot_outfile<-'stat_quant.pdf'  # bar charts
stat_outfile<-'stat_quant.txt'  # tabulated statistical data
pgroups_outfile<-'prot_quant.txt' # filtered set of protein quantitations
drv<-dbDriver('SQLite')
conn<-NULL
dts_names<-NULL
emf<-c('ELF', 'UMTS', 'WiFi') # EM field types
org2cell<-list(mouse = c('IB10'), human = c('U2OS', 'VH10')) # lookup table for cell lines

# function definitions
mybarchart<-function(df, cell_line) {
   print(barchart(paste0(emf, ':', data_set) ~ mad,
      data = df,
      groups = ratio_class,
      #labels = paste0(df$Np, ' / ', df$Nr),
      auto.key = list(space = 'top', columns = 2, title = 'Classes of SILAC ratios',
                      rectangles = TRUE, points = FALSE, cex.title = 1),
      par.settings = simpleTheme(col = c('grey80','grey10')),
      main = paste0('Variability of protein quantitations for ', cell_line, ' cell line'),
      xlab = expression(paste('Variability [MAD(', log[2], '(ratio))]')),
      ylab = 'Data set',
      panel = function(x, y, labels, ...) {
         panel.grid(h = 0, v = -1)
         panel.barchart(x, y, ...)
         #panel.text(x, y, label = labels, cex = 0.6, pos = c(1, 3))
      },
      horizontal = TRUE))
}

mysummarize<-function(df) {
   return(summarise(
      group_by(df, emf, org, cell_line, data_set, ratio_class),
         Np = length(unique(grp_id)), # number of unique protein groups
         Nr = n(),
         min = round(min(log2_ratio), 4),
         max = round(max(log2_ratio), 4),
         median = round(median(log2_ratio), 4),
         mean = round(mean(log2_ratio), 4),
         IQR = round(IQR(log2_ratio), 4),
         sd = round(sd(log2_ratio), 4),
         mad = round(mad(log2_ratio, constant = 1), 4))) # changed the default scale factor (constant)
}

# fetch quantitative data from the database files into a single table
unlink(stat_outfile)
pdf(plot_outfile)

for (o in names(org2cell)) {
   for (c in unlist(org2cell[o])) {
      stats<-NULL
      df<-NULL
      for (e in emf) {
         base_name<-paste0(toupper(e), '_', o)
         dbfile<-file.path(data_dir, base_name, paste0(base_name, '.sqlite'))
         conn<-dbConnect(drv, dbname = dbfile)

         # prepare parts of SQL query used in WHERE clause
         not_null_ratios<-paste0('(', # select protein groups with all six SILAC ratios
            c, '_L0_M0_H1_norm_ratio_HL AND ',
            c, '_L0_M0_H1_norm_ratio_HM AND ',
            c, '_L0_M0_H1_norm_ratio_LM AND ',
            c, '_L1_M1_H0_norm_ratio_LH AND ',
            c, '_L1_M1_H0_norm_ratio_MH AND ',
            c, '_L1_M1_H0_norm_ratio_LM)')

         incons_ratios<-paste0(not_null_ratios, # select protein groups with inconsistent reciprocal ratios
            ' AND NOT ((',
            c, '_L1_M1_H0_norm_ratio_HL > 1 AND ',
            c, '_L1_M1_H0_norm_ratio_HM > 1 AND ',
            c, '_L0_M0_H1_norm_ratio_LH > 1 AND ',
            c, '_L0_M0_H1_norm_ratio_MH > 1) OR (',
            c, '_L1_M1_H0_norm_ratio_LH > 1 AND ',
            c, '_L1_M1_H0_norm_ratio_MH > 1 AND ',
            c, '_L0_M0_H1_norm_ratio_HL > 1 AND ',
            c, '_L0_M0_H1_norm_ratio_HM > 1))')

         cons_ratios<-paste0('((', # select protein groups with consistent reciprocal ratios
            c, '_L1_M1_H0_norm_ratio_HL > 1 AND ',
            c, '_L1_M1_H0_norm_ratio_HM > 1 AND ',
            c, '_L0_M0_H1_norm_ratio_LH > 1 AND ',
            c, '_L0_M0_H1_norm_ratio_MH > 1) OR (',
            c, '_L1_M1_H0_norm_ratio_LH > 1 AND ',
            c, '_L1_M1_H0_norm_ratio_MH > 1 AND ',
            c, '_L0_M0_H1_norm_ratio_HL > 1 AND ',
            c, '_L0_M0_H1_norm_ratio_HM > 1))')

         cons_treat_gt_ctrl_ratios<-paste0('(', # select protein groups for which consistent changes in treatment are
                                                # greater that for control condition
            cons_ratios, ' AND MIN(MAX(',
            c, '_L0_M0_H1_norm_ratio_HL, ',
            c, '_L0_M0_H1_norm_ratio_LH), MAX(',
            c, '_L0_M0_H1_norm_ratio_HM, ',
            c, '_L0_M0_H1_norm_ratio_MH), MAX(',
            c, '_L1_M1_H0_norm_ratio_LH, ',
            c, '_L1_M1_H0_norm_ratio_HL), MAX(',
            c, '_L1_M1_H0_norm_ratio_MH, ',
            c, '_L1_M1_H0_norm_ratio_HM)) > MAX(', 
            c, '_L0_M0_H1_norm_ratio_LM, ',
            c, '_L0_M0_H1_norm_ratio_ML, ',
            c, '_L1_M1_H0_norm_ratio_LM, ',
            c, '_L1_M1_H0_norm_ratio_ML))')

         where_clause<-list(all = not_null_ratios,
                            cons = cons_ratios,
                            incons = incons_ratios,
                            consTC = cons_treat_gt_ctrl_ratios)

         dts_names<-names(where_clause)
         for (d in dts_names) {
            wc<-unlist(where_clause[d], use.names = FALSE)
            # Note: triplex SILAC experiments with inverse metabolic labeling: *_L0_M0_H1 and *_L1_M1_H0 where
            #   '*' refers to one of cell lines: VH10, U2OS or IB10
            #   Ones (or zeros) after the light (L), medium (M) and heavy (H) labels refer to (no) exposure upon EMF
            # 
            # SILAC ratios classified into two categories according to the experimental conditions:
            #   control = pooled control ratios: L0/M0, L1/M1, including their inverse form  (2 x 2)
            #   treated = pooled treatment ratios: H1/L0, H1/M0, L1/H0, M1/H0, including their inverse form (4 x 2)
            #
            # Select all max. quantitated protein groups (i.e., groups having all SILAC ratios) with(out) additional
            # constraint imposed on these ratio: consistent changes up-/down-regulation AND (T)reated ratio > (C)ontrol ratio
            #
            sql<-paste0("
               SELECT
                  grp_id,
                  ratio_class,
                  ratio_type,
                  ratio_value
               FROM
                  (SELECT
                     grp_id,
                     'treated' ratio_class,
                     'H1L0' ratio_type,
                     ", c, "_L0_M0_H1_norm_ratio_HL ratio_value
                  FROM
                     VVV_PGROUP_QUANT
                  UNION
                  SELECT
                     grp_id,
                     'treated' ratio_class,
                     'L0H1' ratio_type,
                     ", c, "_L0_M0_H1_norm_ratio_LH ratio_value
                  FROM
                     VVV_PGROUP_QUANT
                  UNION
                  SELECT
                     grp_id,
                     'treated' ratio_class,
                     'L1H0' ratio_type,
                     ", c, "_L1_M1_H0_norm_ratio_LH ratio_value
                  FROM
                     VVV_PGROUP_QUANT
                  UNION
                  SELECT
                     grp_id,
                     'treated' ratio_class,
                     'H0L1' ratio_type,
                     ", c, "_L1_M1_H0_norm_ratio_HL ratio_value
                  FROM
                     VVV_PGROUP_QUANT
                  UNION
                  SELECT
                     grp_id,
                     'treated' ratio_class,
                     'H1M0' ratio_type,
                     ", c, "_L0_M0_H1_norm_ratio_HM ratio_value
                  FROM
                     VVV_PGROUP_QUANT
                  UNION
                  SELECT
                     grp_id,
                     'treated' ratio_class,
                     'M0H1' ratio_type,
                     ", c, "_L0_M0_H1_norm_ratio_MH ratio_value
                  FROM
                     VVV_PGROUP_QUANT
                  UNION
                  SELECT
                     grp_id,
                     'treated' ratio_class,
                     'M1H0' ratio_type,
                     ", c, "_L1_M1_H0_norm_ratio_MH ratio_value
                  FROM
                     VVV_PGROUP_QUANT
                  UNION
                  SELECT
                     grp_id,
                     'treated' ratio_class,
                     'H0M1' ratio_type,
                     ", c, "_L1_M1_H0_norm_ratio_HM ratio_value
                  FROM
                     VVV_PGROUP_QUANT
                  UNION
                  SELECT
                     grp_id,
                     'control' ratio_class,
                     'L0M0' ratio_type,
                     ", c, "_L0_M0_H1_norm_ratio_LM ratio_value
                  FROM
                     VVV_PGROUP_QUANT
                  UNION
                  SELECT
                     grp_id,
                     'control' ratio_class,
                     'M0L0' ratio_type,
                     ", c, "_L0_M0_H1_norm_ratio_ML ratio_value
                  FROM
                     VVV_PGROUP_QUANT
                  UNION
                  SELECT
                     grp_id,
                     'control' ratio_class,
                     'L1M1' ratio_type,
                     ", c, "_L1_M1_H0_norm_ratio_LM ratio_value
                  FROM
                    VVV_PGROUP_QUANT
                  UNION
                  SELECT
                     grp_id,
                     'control' ratio_class,
                     'M1L1' ratio_type,
                     ", c, "_L1_M1_H0_norm_ratio_ML ratio_value
                  FROM
                     VVV_PGROUP_QUANT) JOIN VVV_PGROUP_QUANT USING(grp_id)
               WHERE ", wc)
            #cat(sql, '\n')
            cat(e, '\t', c, '\t', d, '\n') # print progress info to STDOUT
            res<-dbSendQuery(conn, sql)
            tbl<-fetch(res, n = -1) # fetch rows
            df<-rbind(df, # append rows to data frame
               data.frame(
                  emf = e,                              # EM field type: ELF, UMTS or WiFi
                  org = o,                              # organism: human or mouse
                  cell_line = c,                        # cell line: VH10, U2OS, IB10
                  data_set = d,                         # (un)filtered protein quantitations sets
                  grp_id = tbl$grp_id,                  # protein group ID
                  ratio_class = tbl$ratio_class,        # SILAC ratio class: control or treatment
                  ratio_type = tbl$ratio_type,          # SILAC ratio type: H1L0, L0H1, H0M1 ...
                  log2_ratio = log2(tbl$ratio_value)))  # log2-transformed SILAC protein ratio
            dbClearResult(res) # clean-up
         }
      }

      # write the filtered set (consTC) of protein quantitation onto tab-delimited text file
      write.table(subset(df, data_set == 'consTC'), file = pgroups_outfile, append = TRUE, row.names = FALSE,
         col.names = TRUE, quote = FALSE, sep = '\t')

      # summarize data
      tbl_stats<-mysummarize(df)
      mybarchart(tbl_stats, c)

      # write summary table onto tab-delimited text file
      write.table(tbl_stats, file = stat_outfile, append = TRUE, row.names = FALSE,
         col.names = TRUE, quote = FALSE, sep = '\t')
   }
}

dev.off()
on.exit(dbDisconnect(conn))
