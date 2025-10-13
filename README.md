# Functions to help access annotations from the re-annotated EPICv2 manifest
https://zenodo.org/records/15181885

### load function:
```
library(devtools)
source_url('https://raw.githubusercontent.com/bethan-mallabar-rimmer/EPICv2_manifest/main/functions.R')
```

### get one row per CpG/GENCODEv47 gene annotation:
```
gene_annotations <- expand_annotation(manifest, by='gene')

head(gene_annotations[,1:4]
#Output:
#           IlmnID            Gene       Name AddressA_ID
#1 cg25324105_BC11          ZNF781 cg25324105     1754126
#2 cg25324105_BC11 ENSG00000267552 cg25324105     1754126
#3 cg25383568_TC11           ACTN4 cg25383568    79792482
#4 cg25383568_TC11 ENSG00000298338 cg25383568    79792482
#5 cg25455143_BC11            MBD3 cg25455143    80699190
#6 cg25455143_BC11 ENSG00000267059 cg25455143    80699190
```
### get one row per CpG/GENCODEv47 transcript annotation:
```
transcript_annotations <- expand_annotation(manifest, by='transcript')
head(transcript_annotations[,1:4]
#Output:
#           IlmnID      Transcript       Name AddressA_ID
#1 cg25324105_BC11      ZNF781-202 cg25324105     1754126
#2 cg25324105_BC11 ENST00000586606 cg25324105     1754126
#3 cg25324105_BC11      ZNF781-204 cg25324105     1754126
#4 cg25383568_TC11       ACTN4-201 cg25383568    79792482
#5 cg25383568_TC11       ACTN4-202 cg25383568    79792482
#6 cg25383568_TC11       ACTN4-203 cg25383568    79792482
```

### Pathway analysis - get a list of genes (or transcripts) with significant CpGs in them
```
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
