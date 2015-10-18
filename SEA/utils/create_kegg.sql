--
-- SQL batch script to populate KEGG PATHWAY database in SQLite.
--
-- Author: Arnold Kuzniar
--
-- Version: 1.0
--
--
-- create tables
DROP TABLE IF EXISTS GENE2PROT;
CREATE TABLE GENE2PROT
(
    gene_id  TEXT NOT NULL, -- KEGG gene ID
    prot_acc TEXT NOT NULL  -- UniProtKB accession
);

DROP TABLE IF EXISTS GENE2PATH;
CREATE TABLE GENE2PATH
(
    gene_id TEXT NOT NULL, -- KEGG gene ID
    path_id TEXT NOT NULL  -- KEGG pathway ID
);

DROP TABLE IF EXISTS PATH;
CREATE TABLE PATH
(
    path_id   TEXT NOT NULL, -- KEGG pathway ID
    path_name TEXT           -- KEGG pathway name
);

-- load tab-delimited text files into tables
.sep \t
.import kegg_gene2uniprot.txt GENE2PROT
.import kegg_gene2path.txt GENE2PATH
.import kegg_path.txt PATH

-- create indices on table columns
CREATE INDEX idx_GENE2PROT_gene_id ON GENE2PROT(gene_id);
CREATE INDEX idx_GENE2PROT_prot_acc ON GENE2PROT(prot_acc);
CREATE INDEX idx_GENE2PATH_gene2path ON GENE2PATH(gene_id, path_id);
CREATE INDEX idx_PATH_path_id ON PATH(path_id);
CREATE INDEX idx_PATH_path_name ON PATH(path_name);

-- write human/mouse UniProt-to-KEGG pathway mappings into tab-delimited text files
.header on
.output gene_association.kegg_human
SELECT
   path_id,
   path_name,
   GROUP_CONCAT(prot_acc) AS uniprot_acc
FROM
   PATH INNER JOIN GENE2PATH USING(path_id) INNER JOIN GENE2PROT USING(gene_id)
WHERE
   gene_id LIKE 'hsa%'
GROUP BY path_id;

.output gene_association.kegg_mouse
SELECT
   path_id,
   path_name,
   GROUP_CONCAT(prot_acc) AS uniprot_acc
FROM
   PATH INNER JOIN GENE2PATH USING(path_id) INNER JOIN GENE2PROT USING(gene_id)
WHERE
   gene_id LIKE 'mmu%'
GROUP BY path_id;
