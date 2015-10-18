========================================
Differential protein expression analysis
========================================

Methods:
  fcsigB.py   - fold-change (FC) + peak intensity-based significance B (sigB)
  mzscore.py  - outlier detection based on Z- or M-scores
  fcros.py    - rank-based method (fcros)
  rankprod.py - rank-based method (RankProd)
  limma.R     - linear modeling with empirical Bayes estimation

EMF exposure data sets on human (VH10 and U2OS) and mouse (IB10) cell lines:
  http://piqmie.biotools.nl/download/ELF_human/ELF_human.sqlite
  http://piqmie.biotools.nl/download/ELF_human/ELF_mouse.sqlite
  http://piqmie.biotools.nl/download/UMTS_human/UMTS_human.sqlite
  http://piqmie.biotools.nl/download/UMTS_mouse/UMTS_mouse.sqlite
  http://piqmie.biotools.nl/download/WIFI_human/WIFI_human.sqlite
  http://piqmie.biotools.nl/download/WIFI_mouse/WIFI_mouse.sqlite

Auxiliary scripts used with the Perseus software (v1.3.0.4):
  pre_perseus_sigB.R
  post_perseus_sigB.R
  
  Note: Both R scripts are required for the FC+sigB analysis (see fcsigB.py).
