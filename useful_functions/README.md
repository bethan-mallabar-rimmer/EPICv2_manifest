# Outdated/work in progress - Functions to help access annotations from the re-annotated EPICv2 manifest
Re-annotated manifest available here: https://doi.org/10.5281/zenodo.14933468

### load functions:
```
library(devtools)
source_url('https://raw.githubusercontent.com/bethan-mallabar-rimmer/EPICv2_manifest/main/functions.R')
```

### Using the function - contents of this page:
- Example usage for pathway analysis
- Getting the GENCODEv47 gene and transcript annotations
  - Filtering the annotation to CpGs in the gene body or <1500bp or <200bp upstream of the TSS
- Getting the distance-based promoter and enhancer annotations
  - Filtering the annotation to promoters or enhancers only
- Getting the GeneHancer promoter and enhancer annotations
  - Filtering the annotation to promoters or enhancers only
- Outputting a list or table of filtered genes

## Example use: get a list of genes for pathway analysis

### Get all CpG sites annotated to a gene body (according to GENCODEv47 database)
```
gene_body_cpgs <- expand_annotation(manifest,by='gene') %>% filter_to_genebody()
```

### Get list of all genes with a significant CpG site located in the gene body (according to GENCODEv47 database)
```
my_significant_cpgs <- c('cg25383568_TC11','cg25595446_BC11','cg25908985_BC11','cg25459778_BC11') #...etc. 
gene_body_sig_cpgs <- expand_annotation(manifest,by='gene') %>% filter_to_genebody(sig_cpgs = my_significant_cpgs)
get_annotated_gene_list(gene_body_sig_cpgs)
#Output: a vector c("ACTN4","ENSG00000298338","PRMT1","IHH")

```
### Count the number of significant CpGs per gene
```
get_annotated_gene_table(gene_body_sig_cpgs)
#Output: a table ACTN4 ENSG00000298338             IHH           PRMT1 
                   1         1                      1               1 
```

### Count the number of all CpGs (significant or not) per gene
```
get_annotated_gene_table(gene_body_cpgs)
```

## Getting the GENCODEv47 gene and transcript annotations
### Get one row per CpG/GENCODEv47 gene annotation:
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

## Getting the distance-based promoter and enhancer annotations
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

