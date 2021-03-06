========================================
Differential protein expression analysis
========================================

Methods:
  fcsigB.py   - fold-change (FC) + peak intensity-based significance B (sigB)
  mzscore.py  - outlier detection based on Z- or M-scores
  fcros.py    - rank-based method (fcros)
  rankprod.py - rank-based method (RankProd)
  limma.R     - linear modeling with empirical Bayes estimation

Auxiliary scripts used with the Perseus software (v1.3.0.4) and required for FC+sigB analysis:
  pre_perseus_sigB.R
  post_perseus_sigB.R

EMF exposure data on human (VH10 and U2OS) and mouse (IB10) cell lines are available  for
download from the PIQMIe proteomics server, http://piqmie.biotools.nl/download/...
  ELF_human/ELF_human.sqlite
  ELF_mouse/ELF_mouse.sqlite
  UMTS_human/UMTS_human.sqlite
  UMTS_mouse/UMTS_mouse.sqlite
  WIFI_human/WIFI_human.sqlite
  WIFI_mouse/WIFI_mouse.sqlite

data.tar.bz2 - contains the output files (*.tab) with differentially regulated proteins detected by the different methods
