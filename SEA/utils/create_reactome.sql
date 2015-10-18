--
-- SQL batch script to populate Reactome database in SQLite.
--
-- Author: Arnold Kuzniar
--
-- Version: 1.0
--
--
-- create tables
DROP TABLE IF EXISTS PROT2PATH;
CREATE TABLE PROT2PATH
(
    prot_acc  TEXT NOT NULL, -- UniProtKB accession
    path_id   TEXT NOT NULL, -- Reactome pathway ID
    url       TEXT NOT NULL, -- URL to the Reactome pathway
    path_name TEXT NOT NULL, -- Reactome pathway name
    evidence  TEXT NOT NULL, -- evidence of the association
    org       TEXT NOT NULL  -- organism (species name)
);

-- load tab-delimited text files into tables
.sep \t
.import UniProt2Reactome.txt PROT2PATH

-- create indices on table columns
CREATE INDEX idx_PROT2PATH_prot_acc ON PROT2PATH(prot_acc);
CREATE INDEX idx_PROT2PATH_path_id ON PROT2PATH(path_id);
CREATE INDEX idx_PROT2PATH_org ON PROT2PATH(org);

-- write human/mouse UniProt-to-Reactome pathway mappings into tab-delimited text files
.header on
.output gene_association.reactome_human
SELECT
   path_id,
   path_name,
   GROUP_CONCAT(prot_acc) AS uniprot_acc
FROM
   PROT2PATH
WHERE
   org LIKE 'Homo%'
   -- AND evidence != 'IEA'
GROUP BY path_id;

.output gene_association.reactome_mouse
SELECT
   path_id,
   path_name,
   GROUP_CONCAT(prot_acc) AS uniprot_acc
FROM
   PROT2PATH 
WHERE
   org LIKE 'Mus%'
   -- AND evidence != 'IEA'
GROUP BY path_id;
