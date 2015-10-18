#!/bin/bash
#
# This script retrieves human/mouse mappings from the KEGG pathway database.
#
#
# Author: Arnold Kuzniar
#
# Version: 1.0
#

BASE_URL=http://rest.kegg.jp
ORG=(hsa mmu)
CURL=curl
SFX=txt
GENE2PROTEIN=kegg_gene2uniprot
GENE2PATHWAY=kegg_gene2path
PATHWAY=kegg_path

# delete files from previous run
rm -f $GENE2PROTEIN.$SFX $GENE2PATHWAY.$SFX $PATHWAY.$SFX

# show KEGG pathway release info
$CURL $BASE_URL/info/pathway

for o in ${ORG[@]}; do
   # append file with KEGG geneID-to-UniProtACC mappings
   $CURL $BASE_URL/conv/uniprot/$o | sed 's/up://' >> kegg_gene2uniprot.$SFX
 
   # append file with KEGG geneID-to-pathway  mappings
   $CURL $BASE_URL/link/pathway/$o | sed 's/path://' >> kegg_gene2path.$SFX

   # append file with KEGG pathways
   $CURL $BASE_URL/list/pathway/$o | sed 's/path://' >> kegg_path.$SFX
done

exit 0
