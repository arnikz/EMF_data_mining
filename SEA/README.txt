=======================
Set Enrichment Analysis
=======================

Main script:
  mgsa.R - Bayesian model-based protein set analysis (MGSA)

Databases on pathway/biological process assocations:
  UniProt-GOA release 142 for human; release 128 for mouse (March 2015)
  KEGG PATHWAY release 73.0+/03-10 (March 2015)
  Reactome version 52.0 (March 2015)

Input files:
  {ELF|UMTS|WIFI|RF|EMF}_{human|VH10|U2OS|mouse}_{diff|all}_acc.txt

   diff - study set of differentially regulated proteins (accessions)
   all  - population set of quantitated proteins
   RF   - pooled accessions from UMTS and WIFI 
   EMF  - pooled protein accessions from ELF and RF
   human - pooled accessions from VH10 and U2OS
   mouse - accessions from IB10  

  gene_association.{goa|kegg|reactome}_{human|mouse} - gene/protein association files

Output files:
  mgsa_{goa|kegg|reactome}.txt - table of marginal probability estimates for pathways/biological processes
  mgsa_{goa|kegg|reactome}.pdf - plots of MGSA model fits

Auxiliary scripts:
  fetch_mgsa_pset.py         - fetches a population set required for the MGSA analysis
  fetch_kegg.sh              - fetches human/mouse mappings from the KEGG PATHWAY database 
  create_{kegg|reactome}.sql - SQL batch scripts to populate KEGG or Reactome database locally in SQLite
                             - writes human/mouse pathway association files

  sqlite3 kegg.sqlite < create_kegg.sql
  sqlite3 reactome.sqlite < create_reactome.sql

