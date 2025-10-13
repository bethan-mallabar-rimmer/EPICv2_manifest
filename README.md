# EPICv2_manifest
Functions to help access annotations from the re-annotated EPICv2 manifest: https://zenodo.org/records/15181885

```
library(devtools)
source_url('https://raw.githubusercontent.com/bethan-mallabar-rimmer/EPICv2_manifest/main/functions.R')

#get one row per CpG/GENCODEv47 gene annotation:
gene_annotations <- expand_annotation(manifest, by='gene')

#get one row per CpG/GENCODEv47 transcript annotation:
transcript_annotations <- expand_annotation(manifest, by='transcript')

#pathway analysis - get a list of genes (or transcripts) with significant CpGs in them
my_significant_cpgs <- c('cg25383568_TC11','cg25595446_BC11','cg25908985_BC11') #...etc.
sig_gene_annotations <- gene_annotations[gene_annotations[,1] %in% my_significant_cpgs,]
unique(sig_gene_annotations$Gene)
#output: "ACTN4" "ENSG00000298338" "PRMT1" "IHH"

#get number of significant CpGs per gene
table(sig_gene_annotations$Gene)
#output: ACTN4 ENSG00000298338             IHH           PRMT1
#          1          1                     1               1

#similar for transcripts
sig_transcript_annotations <- transcript_annotations[transcript_annotations[,1] %in% my_significant_cpgs,]
unique(sig_transcript_annotations$Transcript)
#output: "ACTN4-201"       "ACTN4-202"       "ACTN4-203"       "ACTN4-204"       "ACTN4-211"      
#"ENST00000754980" "ENST00000754981" "ENST00000754982" "PRMT1-201"       "PRMT1-202"      
#"PRMT1-217"       "PRMT1-221"       "IHH-201"
table(sig_transcript_annotations$Transcript)
#you get the idea

```
