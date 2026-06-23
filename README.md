## Release Notes

**reannotation.R** is the reannotation code for version 3.0 of the EPICv2 reannotated manifest by Mallabar-Rimmer et al., available from https://doi.org/10.5281/zenodo.14933468 and released on 15th June 2026.

-----------------------------------------------------------------------------------------------------------------------------------------
### Major changes in version 3.0:
- Updated to GENCODEv49
Meaning all column names beginning with 'GENCODEv47...' now begin with 'GENCODEv49' instead. There are no huge differences between these GENCODE releases (same number of genes etc.) but as we were making other updates that required reuploading the manifest on Zenodo, thought we might as well update to the latest version.

- TSS1500 and TSS200 annotations are mutually exclusive
To keep our annotation consistent with Illumina's definition of TSS200 and TSS1500 (https://knowledge.illumina.com/microarray/general/microarray-general-reference_material-list/000009183), the TSS200 annotation covers 1-200bp downstream of the TSS, and the TSS1500 covers 201-1500bp downstream.
In the previous version of our reannotation (2.0), the TSS200 covered 1-200bp downstream and the TSS1500 covered 1-1500bp downstream, so there was overlap.

- The columns beginning with DB (DB_Element_Type", "DB_Element_Gene_Name", "DB_Element_Gene_ID", "DB_Element_Gene_Type") have been removed.
The DB (distance-based) annotation columns in versions 1.0 and 2.0 labelled CpG sites between 500bp downstream and 5000bp upstream of a transcription start site (TSS). The annotation was added by a collaborator based on their previous work in different epigenetic data types, who has requested we remove the annotation.
Our original reannotation used GENCODE version 47 data. I had assumed the "DB" columns also used GENCODEv47, but it was recently clarified that these columns use GENCODEv41. Apologies for the error in describing these columns.

### Other changes
- In the "GENCODEv49_Feature_Type_Specific" column, the wording "exon[number]_CDS" has been changed to "exon[number]_in_CDS", to match the formatting of other annotations in this column for non-coding genes. Version 2.0 used the wording "exon[number]_in_RNA", "exon[number]_in_pseudogene", and "exon[number]_CDS" (so not very consistent).

- In the "GENCODEv49_Feature_Type" column:
Previously introns in the UTR were annotated as "5UTR_intron" and "3UTR_intron", whereas exons in the UTR were annotated as "5UTR" and "3UTR". Exons in UTRs are now annotated as "5UTR_exon" and "3UTR_exon" to make this clearer.

-----------------------------------------------------------------------------------------------------------------------------------------

## Documentation

Run in R programming language, version 4.3.2.

The code has been divided into sections using comments.

### Sections 1-10 - add GENCODEv49 annotation
Input: GENCODEv49 basic gene set from https://www.gencodegenes.org/human/release_49.html and Illumina EPICv2 manifest downloaded from https://support.illumina.com.cn/downloads/infinium-methylationepic-v2-0-product-files.html
- Reformat/clean up
- Use genomic position data from GENCODEv49 to add the following annotations not included in GENCODEv49: introns, UTR directionality (5' and 3'), and TSS200/1500
- Annotate individual CpG sites with info from GENCODE
- Merge annotation back into the original manifest file
Output: manifest file annotated with GENCODEv49 data

### Section 11 - add GeneHancer annotation
Input: above manifest file and GeneHancer data - must be downloaded separately for each chromosome from https://genome-euro.ucsc.edu/cgi-bin/hgTables, see code section 11A 
- Annotate individual CpG sites with info from GeneHancer
- Add true/false column of whether each site is in any GH promoter or enhancer
- Merge this annotation into manifest file and remove GeneHancer data that we cannot redistribute
Output: manifest file annotated with GENCODEv49 & GeneHancer data

### Section 12 - add Horvath and MethylDetectR clock annotation
Input: above manifest file with GENCODEv49 & GeneHancer data, Horvath clock sites downloaded from https://dnamage.clockfoundation.org/importanthints-page  , MethylDetectR clock sites downloaded from https://zenodo.org/records/7154750 
- Add two true/false columns indicating whether sites in manifest are required as input to either clock
Output: fully annotated manifest file requiring cleanup and reformatting

### Sections 13-15 - cleanup/final formatting
Input: above manifest file
- Add back to manifest control probes which did not require re-annotating
- Remove redundant Illumina columns
- Export
Output: our reannotated EPICv2 manifest version 3 available from https://doi.org/10.5281/zenodo.14933468
