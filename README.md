# Functions to help access annotations from the re-annotated EPICv2 manifest
Re-annotated manifest available here: https://zenodo.org/records/15181885

### load function:
```
library(devtools)
source_url('https://raw.githubusercontent.com/bethan-mallabar-rimmer/EPICv2_manifest/main/functions.R')
```

### get one row per CpG/GENCODEv47 gene annotation:
```
gene_annotations <- expand_annotation(manifest, by='gene')

head(gene_annotations[,1:4]
#Output: (many columns excluded)
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
#Output: (many columns excluded)
#           IlmnID      Transcript       Name AddressA_ID
#1 cg25324105_BC11      ZNF781-202 cg25324105     1754126
#2 cg25324105_BC11 ENST00000586606 cg25324105     1754126
#3 cg25324105_BC11      ZNF781-204 cg25324105     1754126
#4 cg25383568_TC11       ACTN4-201 cg25383568    79792482
#5 cg25383568_TC11       ACTN4-202 cg25383568    79792482
#6 cg25383568_TC11       ACTN4-203 cg25383568    79792482
```
### get one row per CpG/distance-based promoter and enhancer annotation:
Promoter = region 1500bp upstream to 500bp downstream of a gene's transcription start site.
Enhancer = region 5000bp upstream of gene's transcription start site, including intergenic regions only (so excluding any region that overlapped the gene body of another gene).
```
regulatory_annotations <- expand_annotation(manifest, by='db')

head(regulatory_annotations[,1:5]
#Output: (many columns excluded)
#           IlmnID   Gene  Element       Name AddressA_ID
#1 cg25324105_BC11 ZNF781 Promoter cg25324105     1754126
#2 cg25383568_TC11                 cg25383568    79792482
#3 cg25455143_BC11                 cg25455143    80699190
#4 cg25459778_BC11 CRAMP1 Promoter cg25459778    60797262
#5 cg25487775_BC11                 cg25487775     5799427
#6 cg25595446_BC11                 cg25595446    65640459
```

### Pathway analysis - get a list of genes with significant CpGs in the gene body, or <1500bp upstream of the TSS:
```
my_significant_cpgs <- c('cg25383568_TC11','cg25595446_BC11','cg25908985_BC11','cg25459778_BC11') #...etc.

sig_gene_annotations <- gene_annotations[gene_annotations[,1] %in% my_significant_cpgs,]

#Get gene list:
unique(sig_gene_annotations$Gene)
#output: "ACTN4" "ENSG00000298338" "CRAMP1" "PRMT1" "IHH" 

#Get number of significant CpGs per gene:
table(sig_gene_annotations$Gene)
#output: ACTN4 ENSG00000298338    CRAMP1       PRMT1      IHH
#          1          1              1           1         1
```
### Get a list of transcripts with significant CpGs in them:
```
sig_transcript_annotations <- transcript_annotations[transcript_annotations[,1] %in% my_significant_cpgs,]

#Get transcript list:
unique(sig_transcript_annotations$Transcript)
#output: [1] "ACTN4-201"       "ACTN4-202"       "ACTN4-203"       "ACTN4-204"       "ACTN4-211"      
#[6] "ENST00000754980" "ENST00000754981" "ENST00000754982" "CRAMP1-201"      "CRAMP1-202"     
#[11] "PRMT1-201"       "PRMT1-202"       "PRMT1-217"       "PRMT1-221"       "IHH-201" 

#Get number of significant CpGs per transcript:
table(sig_transcript_annotations$Transcript)
#you get the idea
```
### Get a list of genes with significant CpGs in their promoter or enhancer (with promoter/enhancer defined according to distance from TSS):
```
sig_regulatory_annotations <- regulatory_annotations[regulatory_annotations[,1] %in% my_significant_cpgs,]

#Get gene list:
unique(sig_regulatory_annotations$Gene)
#output: ""    "CRAMP1"

#Get number of significant CpGs located in a promoter and/or enhancer of each gene:
table(sig_regulatory_annotations$Gene)
#       CRAMP1       - blank space shows 3 out of the 4 CpGs in my_significant_cpgs were not in the promoter or enhancer of any gene
     3      1 
```
