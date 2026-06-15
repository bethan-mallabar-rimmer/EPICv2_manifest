#####################
### RELEASE NOTES ###
#####################

#This is reannotation code for version 3.0 of the EPICv2 reannotated manifest by Mallabar-Rimmer et al., released on [date]

###########################################
#!#!# Major changes in version 3.0: #!#!#
###########################################
# - Updated to GENCODEv49
#   -----------------------------------------------------------------------------------------------------------------------------------------
#   Meaning all column names beginning with 'GENCODEv47...' now begin with 'GENCODEv49' instead.
#   There are no huge differences between these GENCODE releases (same number of genes etc.) but as we were making other updates that required
#   reuploading the manifest on Zenodo, thought we might as well update to the latest version.

# - TSS1500 and TSS200 annotations are mutually exclusive
#   -----------------------------------------------------------------------------------------------------------------------------------------
#   To keep our annotation consistent with Illumina's definition of TSS200 and TSS1500 (https://knowledge.illumina.com/microarray/general/microarray-general-reference_material-list/000009183),
#   the TSS200 annotation covers 1-200bp downstream of the TSS, and the TSS1500 covers 201-1500bp downstream. 
#   In the previous version of our reannotation (2.0), the TSS200 covered 1-200bp downstream and the TSS1500 covered 1-1500bp downstream, so there was overlap.

# - The columns beginning with DB (DB_Element_Type", "DB_Element_Gene_Name", "DB_Element_Gene_ID", "DB_Element_Gene_Type") have been removed.
#   -----------------------------------------------------------------------------------------------------------------------------------------
#   The DB (distance-based) annotation columns in versions 1.0 and 2.0 labelled CpG sites between 500bp downstream and 5000bp upstream of a
#   transcription start site (TSS).
#   The annotation was added by a collaborator based on their previous work in different epigenetic data types, who has requested we remove the
#   annotation.
#   Our original reannotation used GENCODE version 47 data. I had assumed the "DB" columns also used GENCODEv47, but it was recently clarified
#   that these columns use GENCODEv41. Apologies for the error in describing these columns.

### Other changes ###
#- In the "GENCODEv49_Feature_Type_Specific" column, the wording "exon[number]_CDS" has been changed to "exon[number]_in_CDS",
#  to match the formatting of other annotations in this column for non-coding genes. Version 2.0 used the wording
#  "exon[number]_in_RNA", "exon[number]_in_pseudogene", and "exon[number]_CDS" (so not very consistent).

#- In the "GENCODEv49_Feature_Type" column:
    #Previously introns in the UTR were annotated as "5UTR_intron" and "3UTR_intron", whereas exons in the UTR were annotated as "5UTR" and "3UTR".
    #Exons in UTRs are now annotated as "5UTR_exon" and "3UTR_exon" to make this clearer.


#########################
### REANNOTATION CODE ###
#########################
#=========================
#1. Import GENCODEv49 data
#=========================
#GENCODE Human release 49 basic gene annotation, for genome version GRCh38.p14
#Downloaded from the GENCODE website
#In terminal:
#wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz

#workdir <- [working file path not ending in /]

library(rtracklayer)
g <- import(paste0(workdir, "/GENCODE/gencode.v49.basic.annotation.gtf.gz"))

gd <- data.frame(CHR=as.character(seqnames(g)),
                 start=start(g),
                 end=end(g),
                 strand=as.character(strand(g)),
                 gene_name=g$gene_name,
                 gene_id=g$gene_id,
                 transcript_name=g$transcript_name,
                 transcript_id=g$transcript_id,
                 gene_type=g$gene_type,
                 type=as.character(g$type),
                 exon_number=g$exon_number)

#Add an exon number column formatted as e.g. "exon1" rather than 1
gd$exon_number_v2 <- rep(NA, nrow(gd))
gd$exon_number_v2[!is.na(gd$exon_number)] <- paste0('exon',gd$exon_number[!is.na(gd$exon_number)])

#gd$type contains "gene"           "transcript"     "exon"           "CDS"           
#"start_codon"    "stop_codon"     "UTR"            "Selenocysteine"
#Remove Selenocysteine, not really relevant (and overlaps with other features unnecessarily
#complicating annotation)
nrow(gd) #5868512
gd <- gd[!(gd$type=='Selenocysteine'),]
nrow(gd) #5868405

#========================================================================================================================
#2. Annotate intra-exonic features (CDS, UTR, start and stop codons - available in GENCODE for protein coding genes only)
#========================================================================================================================
table(gd$type)
#CDS          exon         gene    start_codon  stop_codon  transcript    UTR 
#2048007     2525461       78691      187041      186872      280000     562333 
gd$exon_feature_type <- rep(NA,nrow(gd))
gd$exon_feature_type[gd$type %in% c('CDS','UTR','start_codon','stop_codon')] <- as.character(gd$type[gd$type %in% c('CDS','UTR','start_codon','stop_codon')])
table(gd$exon_feature_type)
#CDS       start_codon   stop_codon    UTR 
#2048007      187041      186872      562333 

#--------------------------------------------------------------------------------------------------------------------------------
#Non coding genes aren't labelled with intra-exonic features, add annotation to account for this, and group non-coding gene types
#--------------------------------------------------------------------------------------------------------------------------------
protein_coding_gene_types <- c(unique(gd$gene_type[!is.na(gd$exon_feature_type) & !(gd$type %in% c('gene','transcript','exon'))]))
paste(protein_coding_gene_types, collapse=', ') #"protein_coding, IG_V_gene, IG_C_gene, IG_J_gene, TR_C_gene, TR_J_gene, TR_V_gene, TR_D_gene, IG_D_gene"
any(is.na(gd$exon_feature_type[(gd$gene_type %in% protein_coding_gene_types) & !(gd$type %in% c('gene','transcript','exon'))])) #FALSE - good - all protein coding genes have intra-exonic features labelled

noncoding_gene_types <- c(unique(gd$gene_type[is.na(gd$exon_feature_type) & !(gd$gene_type %in% protein_coding_gene_types)]))
paste(noncoding_gene_types, collapse=', ')
#"lncRNA, transcribed_unprocessed_pseudogene, miRNA, processed_pseudogene, transcribed_processed_pseudogene, snRNA, unprocessed_pseudogene, misc_RNA, TEC,
#transcribed_unitary_pseudogene, snoRNA, scaRNA, rRNA_pseudogene, unitary_pseudogene, rRNA, IG_V_pseudogene, sRNA, ribozyme, translated_processed_pseudogene,
#vault_RNA, TR_V_pseudogene, IG_C_pseudogene, TR_J_pseudogene, IG_J_pseudogene, IG_pseudogene, artifact, Mt_tRNA, Mt_rRNA"

gd$exon_feature_type[!(is.na(gd$exon_number)) & 
                       gd$gene_type %in% c('rRNA','misc_RNA',
                                           'lncRNA','scRNA','snoRNA','miRNA',
                                           'ribozyme','snRNA','scaRNA','sRNA',
                                           'Mt_rRNA','Mt_tRNA','vault_RNA')] <- 'RNA_exon' #scRNA was in v47 but not v49, no harm in listing it here anyway
gd$exon_feature_type[!(is.na(gd$exon_number)) & 
                       gd$gene_type %in% c('rRNA_pseudogene',
                                           'transcribed_unprocessed_pseudogene',
                                           'processed_pseudogene',
                                           'transcribed_processed_pseudogene',
                                           'transcribed_unitary_pseudogene',
                                           'unprocessed_pseudogene',
                                           'unitary_pseudogene',
                                           'IG_V_pseudogene',
                                           'IG_pseudogene','IG_C_pseudogene',
                                           'IG_J_pseudogene',
                                           'translated_processed_pseudogene',
                                           'TR_J_pseudogene','TR_V_pseudogene')] <- 'pseudogene_exon'
gd$exon_feature_type[!(is.na(gd$exon_number)) & 
                       gd$gene_type == 'TEC'] <- 'TEC_exon'
gd$exon_feature_type[!(is.na(gd$exon_number)) & 
                       gd$gene_type == 'artifact'] <- 'artifact_exon'

table(gd$exon_feature_type)
#artifact_exon   CDS      pseudogene_exon     RNA_exon       start_codon     stop_codon   TEC_exon    UTR 
#109           2048007       28654              253738          187041         186872       1029    562333

table(gd$type[is.na(gd$exon_feature_type)])
#exon         gene   transcript 
#2241931      78691     280000 

paste(unique(gd$gene_type[is.na(gd$exon_feature_type) & gd$type == 'exon']), collapse=', ')
#"protein_coding, IG_V_gene, IG_C_gene, IG_J_gene, TR_C_gene, TR_J_gene, TR_V_gene, TR_D_gene, IG_D_gene" - all coding, as expected


#============================================================
#3. Calculate position of introns based on position of exons:
#============================================================
library(dplyr)
exons_data <- gd %>% filter(type == 'exon')

#Takes about 5 mins:
introns_list <- exons_data %>%
  group_by(transcript_name) %>%
  arrange(start, .by_group = TRUE) %>%
  do({
    x <- . #x is the df (or rather, tibble) including all exons in transcript, arranged by start site
    if (nrow(x) >=2) {
      data.frame(
        CHR = x$CHR[1],
        start = x$end[1:(nrow(x)-1)] + 1,  # Introns start after the previous exon end
        end = x$start[2:nrow(x)] - 1,      # Introns end before the next exon start
        strand = x$strand[1],
        gene_name = x$gene_name[1],
        gene_id = x$gene_id[1],
        transcript_name = x$transcript_name[1],
        transcript_id= x$transcript_id[1],
        gene_type = x$gene_type[1],
        type = 'intron',
        exon_number = 'intron',
        exon_feature_type = 'intron'
      )
    } else { #if no introns
      data.frame(
        CHR = NA,
        start = NA,
        end = NA,
        strand = NA,
        gene_name = NA,
        gene_id = NA,
        transcript_name = NA,
        transcript_id = NA,
        gene_type = NA,
        type = NA,
        exon_number = NA,
        exon_feature_type = NA
      )
    }
  })
introns_list <- introns_list[!is.na(introns_list$CHR),]

nrow(gd) #5868405
gd <- bind_rows(gd, introns_list)
nrow(gd) #8113866

save(gd, introns_list, file=paste0(workdir,'/annotation_files.RData'))

#==========================================================================
#4. Annotate exons and introns with transcript and gene they are located in
#==========================================================================
#get gene start/end coordinates
gd_genes <- gd[gd$type == 'gene',c('gene_name','gene_id','start','end','CHR')] #gene id required as some genes in different locations have the same names
colnames(gd_genes) <- c('gene_name','gene_id','gene_start','gene_end','gene_CHR')

#get transcript start/end coordinates
gd_transcripts <- gd[gd$type == 'transcript',c('transcript_name','transcript_id','start','end','CHR')]
colnames(gd_transcripts) <- c('transcript_name','transcript_id','transcript_start','transcript_end','transcript_CHR')

#add coordinates of genes and transcripts to each feature in gd
unique(gd$type[!is.na(gd$exon_feature_type)]) #"exon" "CDS" "start_codon" "stop_codon" "UTR" "intron" 

gd2 <- gd[!is.na(gd$exon_feature_type),] %>%
  left_join(gd_transcripts[,colnames(gd_transcripts) != 'transcript_id'], by = "transcript_name") %>%
  left_join(gd_genes[,colnames(gd_genes) != 'gene_name'], by = "gene_id")
gd2 <- gd2 %>%
  group_by(transcript_id) %>%
  arrange(start, end, .by_group=TRUE) #NOTE TO SELF - CHANGED transcript_name TO transcript_id - CHECK FOR DOWNSTREAM EFFECTS!

#======================================
#5. Calculate whether UTRs are 5' or 3' - UPDATE on 02/04/2025, included in v2 of reannotation
#======================================
#5' UTRs = upstream, closest to transcription start site (TSS)
#3' UTRs = downstream, closest to transcription end site (TES)

gd2_minus <- gd2[gd2$strand == '-',] #get features on minus strand (TES <- 3'UTR <- introns/exons <- 5'UTR <- TSS)
gd2_plus <- gd2[gd2$strand == '+',] #get features on plus strand (TSS -> 5'UTR -> introns/exons -> 3'UTR -> TES)

#Get the feature nearest to the TSS for each transcript:
gd2_minus_test <- gd2_minus %>%
  group_by(transcript_id) %>%
  arrange(start, end) %>%
  summarise(last(exon_feature_type), .groups = 'keep')
unique(gd2_minus_test[,2]) #Features nearest TSS include misc exons, UTRs, start_codons
table(gd2_minus_test[,2])
#artifact_exon     CDS      pseudogene_exon    RNA_exon     start_codon    TEC_exon        UTR 
#   10             168            7710           37871         1192         454           90288 

gd2_plus_test <- gd2_plus %>%
  group_by(transcript_id) %>%
  arrange(start, end) %>%
  summarise(first(exon_feature_type), .groups = 'keep')
unique(gd2_plus_test[,2]) #Also misc exons, UTRs, start_codons
table(gd2_plus_test[,2])
#artifact_exon   CDS  pseudogene_exon    RNA_exon    start_codon  TEC_exon   UTR 
#     9          189      7498            38761        1373         565      93912

#Strange that some transcripts start with exon rather than UTR, but at least none start on a stop codon.

#Do all transcripts have 2 UTRs?
gd2_utr_count <- gd2 %>%
  group_by(transcript_id) %>%
  summarise(UTR_count = sum(exon_feature_type == 'UTR'))
any(duplicated(gd2_utr_count$transcript_id)) #FALSE
signif(sum(gd2_utr_count$UTR_count == 0)/nrow(gd2_utr_count)*100, 4) #33.25% of transcripts have no UTRs
signif(sum(gd2_utr_count$UTR_count < 2)/nrow(gd2_utr_count)*100, 4) #34.2% have less than the 2 UTRs one at each end which you might expect (this was way higher in v47 - 60%!)
signif(sum(gd2_utr_count$UTR_count > 2)/nrow(gd2_utr_count)*100, 4) #39.68 have >2 UTRs per unique transcript
#Looking at example transcript ENST00000000412.8 (3 UTRs)
as.data.frame(gd2[gd2$transcript_id == 'ENST00000000412.8',c('transcript_id','type','exon_number_v2')])
#transcript_id             type       exon_number_v2
#1  ENST00000000412.8         UTR          exon7
#2  ENST00000000412.8  stop_codon          exon7
#3  ENST00000000412.8         CDS          exon7
#4  ENST00000000412.8      intron           <NA>
#5  ENST00000000412.8         CDS          exon6
#6  ENST00000000412.8      intron           <NA>
# ...
#13 ENST00000000412.8         CDS          exon2
#14 ENST00000000412.8 start_codon          exon2
#15 ENST00000000412.8         UTR          exon2
#16 ENST00000000412.8      intron           <NA>
#17 ENST00000000412.8         UTR          exon1
#This shows the reason for some transcripts having >2 UTRs labelled is because some UTRs span more than one exon (e.g. the UTR above spanning exons 1 and 2)
#I.e. in this case, it looks like exons 1 and 2 are transcribed, but exon 1 and the first part of exon 2 are not translated.

#This could be recorded as e.g.:
#UTR (exon 1), intron, UTR (exon 2), intron, CDS (exon 3), intron, CDS (exon 4), intron, CDS (exon 5), ... UTR (exon 7) 

#UTR count of protein coding vs noncoding genes:
protein_coding_transcripts <- unique(gd2$transcript_id[gd2$gene_type %in% protein_coding_gene_types])
#coding:
sum(gd2_utr_count[gd2_utr_count$transcript_id %in% protein_coding_transcripts,]$UTR_count == 0) #234
sum(gd2_utr_count[gd2_utr_count$transcript_id %in% protein_coding_transcripts,]$UTR_count > 0) #186888
#noncoding:
sum(gd2_utr_count[!(gd2_utr_count$transcript_id %in% protein_coding_transcripts),]$UTR_count == 0) #92878
sum(gd2_utr_count[!(gd2_utr_count$transcript_id %in% protein_coding_transcripts),]$UTR_count > 0) #0
#This shows GENCODE does not annotate UTRs in noncoding genes (makes sense - there can't be an untranslated region if the gene is not translated at all)
#and this accounts for the substantial number of transcripts ending on an exon
#There are also a few exceptions of protein coding transcripts not having UTRs.

#In conclusion:
#the definition of a 5' UTR is that it's between the TSS and the first CDS (coding exon)
#the definition of a 3' UTR is that it's between the last CDS and the TES

#Get the first and last CDS for each transcript:
gd2_first_CDS_p <- gd2_plus[gd2_plus$exon_feature_type == 'CDS',] %>%
  group_by(transcript_id) %>%
  arrange(start, end, .by_group = TRUE) %>%
  slice_head(n = 1)
gd2_first_CDS_p <- gd2_first_CDS_p %>% dplyr::rename(first_CDS_start = start,
                                                       first_CDS_end = end)
gd2_last_CDS_p <- gd2_plus[gd2_plus$exon_feature_type == 'CDS',] %>%
  group_by(transcript_id) %>%
  arrange(start, end, .by_group = TRUE) %>%
  slice_tail(n = 1)
gd2_last_CDS_p <- gd2_last_CDS_p %>% dplyr::rename(last_CDS_start = start,
                                                     last_CDS_end = end)
gd2_first_CDS_m <- gd2_minus[gd2_minus$exon_feature_type == 'CDS',] %>%
  group_by(transcript_id) %>%
  arrange(start, end, .by_group = TRUE) %>%
  slice_tail(n = 1)
gd2_first_CDS_m <- gd2_first_CDS_m %>% dplyr::rename(first_CDS_start = start,
                                                       first_CDS_end = end)
gd2_last_CDS_m <- gd2_minus[gd2_minus$exon_feature_type == 'CDS',] %>%
  group_by(transcript_id) %>%
  arrange(start, end, .by_group = TRUE) %>%
  slice_head(n = 1)
gd2_last_CDS_m <- gd2_last_CDS_m %>% dplyr::rename(last_CDS_start = start,
                                                     last_CDS_end = end)

#Use positions of first and last CDS to calculate which UTRs are 5' or 3':
gd2_plus <- left_join(gd2_plus, gd2_first_CDS_p[,c('transcript_id','first_CDS_start','first_CDS_end')])
gd2_plus <- left_join(gd2_plus, gd2_last_CDS_p[,c('transcript_id','last_CDS_start','last_CDS_end')])

gd2_minus <- left_join(gd2_minus, gd2_first_CDS_m[,c('transcript_id','first_CDS_start','first_CDS_end')])
gd2_minus <- left_join(gd2_minus, gd2_last_CDS_m[,c('transcript_id','last_CDS_start','last_CDS_end')])

save(gd, gd2, gd2_plus, gd2_minus, introns_list, file=paste0(workdir,'/annotation_files.RData'))

gd2_plus$UTR <- rep('', nrow(gd2_plus))
gd2_plus$UTR[(gd2_plus$start >= gd2_plus$transcript_start) & 
               (gd2_plus$end < gd2_plus$first_CDS_start)] <- '5UTR'
gd2_plus$UTR[(gd2_plus$start > gd2_plus$last_CDS_end) & 
               (gd2_plus$end <= gd2_plus$transcript_end)] <- '3UTR'
#TSS ('transcript start' column) -> 5UTR -> first CDS -> ... -> last CDS -> 3UTR -> TES ('transcript end' column)

gd2_minus$UTR <- rep('', nrow(gd2_minus))
gd2_minus$UTR[(gd2_minus$start > gd2_minus$first_CDS_end) & 
               (gd2_minus$end <= gd2_minus$transcript_end)] <- '5UTR'
gd2_minus$UTR[(gd2_minus$start >= gd2_minus$transcript_start) & 
                (gd2_minus$end < gd2_minus$last_CDS_start)] <- '3UTR'
#TES ('transcript start' column) <- 3UTR <- last <- ... <- first <- 5UTR <- TSS ('transcript end' column)

nrow(gd2_minus[gd2_minus$UTR=='' & gd2_minus$exon_feature_type == 'UTR',]) #0, all UTRs classified
nrow(gd2_plus[gd2_plus$UTR=='' & gd2_plus$exon_feature_type == 'UTR',]) #0, all UTRs classified

paste(unique(gd2_plus$exon_feature_type), collapse=', ') #"UTR, start_codon, CDS, intron, stop_codon, pseudogene_exon, RNA_exon, TEC_exon, artifact_exon"
paste(unique(gd2_minus$exon_feature_type), collapse=', ')#"UTR, stop_codon, CDS, intron, start_codon, RNA_exon, pseudogene_exon, artifact_exon, TEC_exon"
gd2_plus$exon_feature_type[gd2_plus$exon_feature_type == 'UTR'] <- gd2_plus$UTR[gd2_plus$exon_feature_type == 'UTR']
gd2_plus$exon_feature_type[gd2_plus$exon_feature_type == 'stop_codon'] <- paste0(gd2_plus$UTR[gd2_plus$exon_feature_type == 'stop_codon'], '_stop_codon')
gd2_plus$exon_feature_type[gd2_plus$exon_feature_type == 'intron' & gd2_plus$UTR!=''] <- paste0(gd2_plus$UTR[gd2_plus$exon_feature_type == 'intron' & gd2_plus$UTR!=''], '_intron')
gd2_minus$exon_feature_type[gd2_minus$exon_feature_type == 'UTR'] <- gd2_minus$UTR[gd2_minus$exon_feature_type == 'UTR']
gd2_minus$exon_feature_type[gd2_minus$exon_feature_type == 'stop_codon'] <- paste0(gd2_minus$UTR[gd2_minus$exon_feature_type == 'stop_codon'], '_stop_codon')
gd2_minus$exon_feature_type[gd2_minus$exon_feature_type == 'intron' & gd2_minus$UTR!=''] <- paste0(gd2_minus$UTR[gd2_minus$exon_feature_type == 'intron' & gd2_minus$UTR!=''], '_intron')
paste(unique(gd2_plus$exon_feature_type), collapse=', ') #"5UTR, start_codon, CDS, intron, 3UTR_stop_codon, 3UTR, 5UTR_intron, 3UTR_intron, pseudogene_exon, RNA_exon, TEC_exon, artifact_exon"
paste(unique(gd2_minus$exon_feature_type), collapse=', ')#"3UTR, 3UTR_stop_codon, CDS, intron, start_codon, 5UTR, 5UTR_intron, 3UTR_intron, RNA_exon, pseudogene_exon, artifact_exon, TEC_exon"

#UPDATE on 17/03/2026, included in v3 of reannotation - changing 5UTR and 3UTR to 5UTR_exon and 3UTR_exon
#because "5UTR_exon" and "5UTR_intron" is clearer than "5UTR" and "5UTR_intron"
gd2_plus$exon_feature_type[gd2_plus$exon_feature_type == '5UTR'] <- '5UTR_exon'
gd2_plus$exon_feature_type[gd2_plus$exon_feature_type == '3UTR'] <- '3UTR_exon'
gd2_minus$exon_feature_type[gd2_minus$exon_feature_type == '5UTR'] <- '5UTR_exon'
gd2_minus$exon_feature_type[gd2_minus$exon_feature_type == '3UTR'] <- '3UTR_exon'
paste(unique(gd2_plus$exon_feature_type), collapse=', ') #"5UTR_exon, start_codon, CDS, intron, 3UTR_stop_codon, 3UTR_exon, 5UTR_intron, 3UTR_intron, pseudogene_exon, RNA_exon, TEC_exon, artifact_exon"
paste(unique(gd2_minus$exon_feature_type), collapse=', ')#"3UTR_exon, 3UTR_stop_codon, CDS, intron, start_codon, 5UTR_exon, 5UTR_intron, 3UTR_intron, RNA_exon, pseudogene_exon, artifact_exon, TEC_exon"

save(gd,gd2, gd2_plus, gd2_minus, introns_list, file=paste0(workdir,'/annotation_files.RData'))

#======================================
#5. Add TSS200 as well as TSS1500 - UPDATE on 04/04/2025, included in v2 of reannotation
#======================================
#ANOTHER UPDATE on 17/03/26 included in v3 of reannotation: Illumina define the TSS1500
#as between 200-1500bp upstream of a TSS (https://knowledge.illumina.com/microarray/general/microarray-general-reference_material-list/000009183),
#whereas in v2 of reannotation we previously defined it as 1-1500bp upstream.
#The new definition in v3 of reannotation is:
#TSS200 = 1-200bp upstream
#TSS1500 = 201-1500bp upstream

gd2_plus_tss200 <- gd2_plus %>%
  group_by(transcript_id) %>%
  arrange(start, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  mutate(start = transcript_start-201, #Don't know the exact details of how Illumina made their TSS200 annotation, but the region between N-1 base pairs and N-201 base pairs encompasses 200 basepairs, which seems sensible for TSS200
         end = transcript_start-1,
         type = 'TSS200',
         exon_feature_type = 'TSS200',
         exon_number = 'TSS200',
         exon_number_v2 = 'TSS200')

gd2_minus_tss200 <- gd2_minus %>%
  group_by(transcript_id) %>%
  arrange(start, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  mutate(start = transcript_end+1,
         end = transcript_end+201,
         type = 'TSS200',
         exon_feature_type = 'TSS200',
         exon_number = 'TSS200',
         exon_number_v2 = 'TSS200')
  
gd2_plus_tss1500 <- gd2_plus %>%
  group_by(transcript_id) %>%
  arrange(start, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  mutate(start = transcript_start-1501, #Start at -202 so there is no overlap with the TSS200
         end = transcript_start-202,
         type = 'TSS1500',
         exon_feature_type = 'TSS1500',
         exon_number = 'TSS1500',
         exon_number_v2 = 'TSS1500')

gd2_minus_tss1500 <- gd2_minus %>%
  group_by(transcript_id) %>%
  arrange(start, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  mutate(start = transcript_end+202,
         end = transcript_end+1501,
         type = 'TSS1500',
         exon_feature_type = 'TSS1500',
         exon_number = 'TSS1500',
         exon_number_v2 = 'TSS1500')

nrow(gd2_minus) + nrow(gd2_plus) == nrow(gd2) #TRUE, good, nothing has been omitted
gd3 <- bind_rows(gd2_plus, gd2_plus_tss200, gd2_plus_tss1500,
                 gd2_minus, gd2_minus_tss200, gd2_minus_tss1500)

#turns out all stop codons are in the 3UTR so no need to label them as such
#users who want to look at the 3'UTR can decide if they want to include stop codons or not, e.g.:
#grep 'stop_codon' and '3UTR' to include them
#grep just '3UTR' to exclude them
gd3$exon_feature_type[gd3$exon_feature_type == '3UTR_stop_codon'] <- 'stop_codon'

save(gd, gd2, gd3, file=paste0(workdir,'/annotation_files.RData'))

#========================================================================================
#6. Ensure there is only one annotation per feature (i.e. no exactly overlapping features
#with the same start and end positions annotated to the same transcript)
#========================================================================================

#Takes like 20 mins
gd3 <- gd3 %>%
  group_by(transcript_id, start, end) %>%
  summarise(
    GENCODEv49_Gene_Name = paste(unique(gene_name), collapse = '; '),
    GENCODEv49_Gene_ID = paste(unique(gene_id), collapse = '; '),
    GENCODEv49_Gene_CHR = paste(unique(CHR), collapse = '; '),
    GENCODEv49_Gene_Strand = paste(unique(strand), collapse = '; '),
    GENCODEv49_Gene_Start = paste(unique(gene_start), collapse='; '),
    GENCODEv49_Gene_End = paste(unique(gene_end), collapse='; '),
    GENCODEv49_Gene_Type = paste(unique(gene_type), collapse = '; '),
    GENCODEv49_Transcript_Name = paste(unique(transcript_name), collapse = '; '),
    GENCODEv49_Transcript_ID = paste(unique(transcript_id), collapse = '; '),
    GENCODEv49_Transcript_Start = paste(unique(transcript_start), collapse='; '),
    GENCODEv49_Transcript_End = paste(unique(transcript_end), collapse='; '),
    GENCODEv49_Feature_Type_Specific = ifelse(('3UTR_exon' %in% exon_feature_type) & 
                                             ('stop_codon' %in% exon_feature_type), #if a UTR and stop codon cover the exact same 3 bases I'm pretty sure that's just a regular old stop codon, so label as such
                                           paste(unique(paste0(exon_number_v2, '_stop_codon')),collapse=';'), 
                                           ifelse(('CDS' %in% exon_feature_type) & 
                                                    ('start_codon' %in% exon_feature_type), #same with a 3-base-CDS that exactly overlaps a start codon. also known as a start codon.
                                                  paste(unique(paste0(exon_number_v2,'_start_codon')),collapse=';'),
                                                  paste(unique(paste0(exon_number_v2,'_',exon_feature_type)),collapse = ';'))),
    
    GENCODEv49_Feature_Type = ifelse(('3UTR_exon' %in% exon_feature_type) & 
                                    ('stop_codon' %in% exon_feature_type),
                                  'stop_codon',
                                  ifelse(('CDS' %in% exon_feature_type) & 
                                           ('start_codon' %in% exon_feature_type),
                                         'start_codon',
                                         paste(unique(exon_feature_type),collapse = ';'))),
    
    GENCODEv49_Feature_Exon_Number = paste(unique(exon_number_v2), collapse = '; '),
    GENCODEv49_Feature_Start = paste(unique(start), collapse = '; '),
    GENCODEv49_Feature_End = paste(unique(end), collapse = '; '),
    .groups = "keep"
  )
gd3$GENCODEv49_Feature_Exon_Number[gd3$GENCODEv49_Feature_Type=='intron'] <- 'intron'
gd3$GENCODEv49_Feature_Type_Specific[gd3$GENCODEv49_Feature_Type=='intron'] <- 'intron'

#sanity check that each unique combination of transcript name, feature start site, and feature end site,
#contains only one feature (i.e. there are no exactly overlapping features in the same transcript)
sum(grepl(';',gd3$GENCODEv49_Feature_Exon_Number)) #=0, good
sum(grepl(';',gd3$GENCODEv49_Feature_Type)) #=0
sum(grepl(';',gd3$GENCODEv49_Feature_Type_Specific)) #=0

save(gd,gd2,gd3,file=paste0(workdir,'/annotation_files.RData'))

#==============================
#7. Clean up feature annotation
#==============================
#remove decimal points from gene IDs
gd4 <- gd3 %>%
  dplyr::mutate(GENCODEv49_Gene_ID = gsub('\\..*','',GENCODEv49_Gene_ID),
         GENCODEv49_Transcript_ID = gsub('\\..*','',GENCODEv49_Transcript_ID))
#clean up some stuff in the feature type specific column:
gd4$GENCODEv49_Feature_Type_Specific <- gsub('_CDS','_in_CDS',gd4$GENCODEv49_Feature_Type_Specific)
gd4$GENCODEv49_Feature_Type_Specific <- gsub('_pseudogene_exon','_in_pseudogene',gd4$GENCODEv49_Feature_Type_Specific)
gd4$GENCODEv49_Feature_Type_Specific <- gsub('_RNA_exon','_in_RNA',gd4$GENCODEv49_Feature_Type_Specific)
gd4$GENCODEv49_Feature_Type_Specific <- gsub('_TEC_exon','_in_TEC',gd4$GENCODEv49_Feature_Type_Specific)
gd4$GENCODEv49_Feature_Type_Specific <- gsub('_artifact_exon','_in_artifact',gd4$GENCODEv49_Feature_Type_Specific)
gd4$GENCODEv49_Feature_Type_Specific <- gsub('_intron','intron',gd4$GENCODEv49_Feature_Type_Specific)
gd4$GENCODEv49_Feature_Type_Specific <- gsub('_3UTR_exon','_in_3UTR',gd4$GENCODEv49_Feature_Type_Specific)
gd4$GENCODEv49_Feature_Type_Specific <- gsub('_5UTR_exon','_in_5UTR',gd4$GENCODEv49_Feature_Type_Specific)
gd4$GENCODEv49_Feature_Type_Specific <- gsub('NA_5UTRintron','intron_in_5UTR',gd4$GENCODEv49_Feature_Type_Specific)
gd4$GENCODEv49_Feature_Type_Specific <- gsub('NA_3UTRintron','intron_in_3UTR',gd4$GENCODEv49_Feature_Type_Specific)
gd4$GENCODEv49_Feature_Type_Specific[gd4$GENCODEv49_Feature_Type_Specific == 'TSS1500_TSS1500'] <-'TSS1500'
gd4$GENCODEv49_Feature_Type_Specific[gd4$GENCODEv49_Feature_Type_Specific == 'TSS200_TSS200'] <-'TSS200'

save(gd,gd2,gd3,gd4,file=paste0(workdir,'/annotation_files.RData'))

#====================================================================================
#8. Annotate individual CpG sites with genes/transcripts/features they are located in
#====================================================================================
#Import Illumina manifest as data.table - downloaded from https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html Infinium MethylationEPIC v2.0 Product Files
library(data.table)
manifest <- fread(paste0(workdir,'/manifest_EPIC-8v2-0_A1.csv'), skip=7, fill=TRUE, data.table=F)
setDT(manifest)
#Extract ID, chromosome, and position of all sites in manifest
m2 <- manifest[,.(IlmnID,CHR,MAPINFO)]
#Add empty columns to m2 with same column names as gd4
gd4$CHR <- gd4$GENCODEv49_Gene_CHR
setDT(gd4)
colz <- colnames(gd4)[!(colnames(gd4) %in% c('transcript_id','start','end','CHR'))]
m2[, (colz):=character()]

#Annotating CpG sites with genes/transcripts/features based on chromosome and
#genomic position is very computationally intensive.
#Therefore this step was done one chromosome at a time, so it is only necessary
#to match based on position, not chromosome.

#Split our annotation by chromosome:
gdc1 <- gd4[CHR == 'chr1',]
gdc2 <- gd4[CHR == 'chr2',]
gdc3 <- gd4[CHR == 'chr3',]
gdc4 <- gd4[CHR == 'chr4',]
gdc5 <- gd4[CHR == 'chr5',]
gdc6 <- gd4[CHR == 'chr6',]
gdc7 <- gd4[CHR == 'chr7',]
gdc8 <- gd4[CHR == 'chr8',]
gdc9 <- gd4[CHR == 'chr9',]
gdc10 <- gd4[CHR == 'chr10',]
gdc11 <- gd4[CHR == 'chr11',]
gdc12 <- gd4[CHR == 'chr12',]
gdc13 <- gd4[CHR == 'chr13',]
gdc14 <- gd4[CHR == 'chr14',]
gdc15 <- gd4[CHR == 'chr15',]
gdc16 <- gd4[CHR == 'chr16',]
gdc17 <- gd4[CHR == 'chr17',]
gdc18 <- gd4[CHR == 'chr18',]
gdc19 <- gd4[CHR == 'chr19',]
gdc20 <- gd4[CHR == 'chr20',]
gdc21 <- gd4[CHR == 'chr21',]
gdc22 <- gd4[CHR == 'chr22',]
gdcX <- gd4[CHR == 'chrX',]
gdcY <- gd4[CHR == 'chrY',]
gdcM <- gd4[CHR == 'chrM',]

#Split CpGs by chromosome:
mf1 <- m2[CHR == 'chr1',]
mf2 <- m2[CHR == 'chr2',]
mf3 <- m2[CHR == 'chr3',]
mf4 <- m2[CHR == 'chr4',]
mf5 <- m2[CHR == 'chr5',]
mf6 <- m2[CHR == 'chr6',]
mf7 <- m2[CHR == 'chr7',]
mf8 <- m2[CHR == 'chr8',]
mf9 <- m2[CHR == 'chr9',]
mf10 <- m2[CHR == 'chr10',]
mf11 <- m2[CHR == 'chr11',]
mf12 <- m2[CHR == 'chr12',]
mf13 <- m2[CHR == 'chr13',]
mf14 <- m2[CHR == 'chr14',]
mf15 <- m2[CHR == 'chr15',]
mf16 <- m2[CHR == 'chr16',]
mf17 <- m2[CHR == 'chr17',]
mf18 <- m2[CHR == 'chr18',]
mf19 <- m2[CHR == 'chr19',]
mf20 <- m2[CHR == 'chr20',]
mf21 <- m2[CHR == 'chr21',]
mf22 <- m2[CHR == 'chr22',]
mfX <- m2[CHR == 'chrX',]
mfY <- m2[CHR == 'chrY',]
mfM <- m2[CHR == 'chrM',]

#==== Annotate CpGs with genes/transcripts/features one chromosome at a time: =====
#The following will take about 2 hours, if you just leave it running it saves after each chromosome
#so if it crashes it can just be picked up from the last chromosome
mf1_2 <- mf1 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc1[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf2_2 <- mf2 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc2[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf3_2 <- mf3 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc3[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf4_2 <- mf4 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc4[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf5_2 <- mf5 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc5[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf6_2 <- mf6 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc6[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf7_2 <- mf7 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc7[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf8_2 <- mf8 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc8[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf9_2 <- mf9 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc9[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf10_2 <- mf10 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc10[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf11_2 <- mf11 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc11[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf12_2 <- mf12 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc12[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf13_2 <- mf13 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc13[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf14_2 <- mf14 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc14[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf15_2 <- mf15 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc15[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf16_2 <- mf16 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc16[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,
     mf16_2, file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf17_2 <- mf17 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc17[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,
     mf16_2, mf17_2, file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf18_2 <- mf18 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc18[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,
     mf16_2, mf17_2, mf18_2, file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf19_2 <- mf19 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc19[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,
     mf16_2, mf17_2, mf18_2, mf19_2, file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf20_2 <- mf20 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc20[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,
     mf16_2, mf17_2, mf18_2, mf19_2, mf20_2, file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf21_2 <- mf21 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc21[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,
     mf16_2, mf17_2, mf18_2, mf19_2, mf20_2, mf21_2, file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mf22_2 <- mf22 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdc22[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,
     mf16_2, mf17_2, mf18_2, mf19_2, mf20_2, mf21_2,
     mf22_2, file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mfX_2 <- mfX %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdcX[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,
     mf16_2, mf17_2, mf18_2, mf19_2, mf20_2, mf21_2,
     mf22_2, mfX_2, file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mfY_2 <- mfY %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdcY[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,
     mf16_2, mf17_2, mf18_2, mf19_2, mf20_2, mf21_2,
     mf22_2, mfX_2, mfY_2, file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))
mfM_2 <- mfM %>%
  group_by(MAPINFO) %>%
  do({
    summarise(gdcM[start <= .$MAPINFO & end >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,
     mf16_2, mf17_2, mf18_2, mf19_2, mf20_2, mf21_2,
     mf22_2, mfX_2, mfY_2, mfM_2, file=paste0(workdir,'/GENCODE_manifest_by_chr.RData'))

#====== Rejoin into one table: ======
mf_list <- list(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,mf16_2,mf17_2,mf18_2,mf19_2,mf20_2,mf21_2,mf22_2,mfX_2,mfY_2,mfM_2)
chrs <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12',
          'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX',
          'chrY','chrM')
for (i in 1:25) {
  mf_list[[i]]$CHR <- chrs[i]
}
m2 <- bind_rows(mf_list$mf1_2,mf_list$mf2_2,mf_list$mf3_2,mf_list$mf4_2,mf_list$mf5_2,mf_list$mf6_2,mf_list$mf7_2,mf_list$mf8_2,mf_list$mf9_2,mf_list$mf10_2,mf_list$mf11_2,mf_list$mf12_2,mf_list$mf13_2,mf_list$mf14_2,mf_list$mf15_2,mf_list$mf16_2,mf_list$mf17_2,mf_list$mf18_2,mf_list$mf19_2,mf_list$mf20_2,mf_list$mf21_2,mf_list$mf22_2,mf_list$mfX_2,mf_list$mfY_2,mf_list$mfM_2)
m2 <- mf_list[[1]]
for (i in 2:25) {
  m2 <- bind_rows(m2,mf_list[[i]])
}

save(m2,file=paste0(workdir,'/1_GENCODE_manifest_by_chr_combined.RData'))

#=========================================
#9. Add in gene body true/false annotation
#=========================================
m2$GENCODEv49_In_Gene_Body <- ifelse(grepl("exon|intron|UTR|codon", m2$GENCODEv49_Feature_Type), TRUE, FALSE)

save(m2,file=paste0(workdir,'/1_GENCODE_manifest_by_chr_combined.RData'))

#==================================================
#10. Append our CpG annotation to Illumina manifest
#==================================================
#Combine chromosome and genomic position (MAPINFO)
manifest$CM <- paste0(manifest$CHR,'_',manifest$MAPINFO)
m2$CM <- paste0(m2$CHR,'_',m2$MAPINFO)

m2 <- dplyr::left_join(manifest, m2[,c('CM','GENCODEv49_Gene_Name','GENCODEv49_Gene_ID','GENCODEv49_Gene_Type',
                                'GENCODEv49_Transcript_Name','GENCODEv49_Transcript_ID',
                                'GENCODEv49_Feature_Type','GENCODEv49_Feature_Exon_Number',
                                'GENCODEv49_Feature_Type_Specific', 'GENCODEv49_In_Gene_Body',
                                'GENCODEv49_Gene_Start','GENCODEv49_Gene_End','GENCODEv49_Gene_Strand',
                                'GENCODEv49_Transcript_Start','GENCODEv49_Transcript_End',
                                'GENCODEv49_Feature_Start','GENCODEv49_Feature_End')],
                by='CM', keep=FALSE)

#Change/correct some formatting
m2$GENCODEv49_Feature_Exon_Number <- gsub('NA','intron',m2$GENCODEv49_Feature_Exon_Number)
m2$GENCODEv49_Feature_Type <- gsub('CDS','CDS_exon',m2$GENCODEv49_Feature_Type)

save(m2,file=paste0(workdir,'/1_GENCODE_manifest_by_chr_combined.RData'))

#===================================================================================================
####################################################################################################
#===================================================================================================
#11. Annotate whether each CpG site is in any double-elite GeneHancer element (enhancer or promoter)
#===================================================================================================
####################################################################################################
#===================================================================================================
#GeneHancer data was used to annotate each CpG with a binary TRUE/FALSE indicator
#of whether the CpG is in any double-elite GeneHancer element (enhancer or promoter).

#It was necessary to download GeneHancer data for this purpose, but due to restrictions on data
#sharing, this data was not included in the final reannotated manifest.

#=============================
#11A. Download GeneHancer data
#=============================
#It is not possible to download all GeneHancer data at once, however, it can be downloaded from
#UCSC table browser one chromosome at a time, so this step take some manual effort.

#Get start and end coordinates of each chromosome:
library(GenomeInfoDb)
chr_lengths <- seqlengths(Seqinfo(genome="hg38"))
chr_lengths <- chr_lengths[paste0("chr", c(1:22, "X", "Y"))]
paste0(names(chr_lengths), ":1-", format(chr_lengths, big.mark=","))

#Output:
#chr1:1-248,956,422
#chr2:1-242,193,529
#chr3:1-198,295,559
#...
#etc.

#GeneHancer data for each chromosome was downloaded from https://genome-euro.ucsc.edu/cgi-bin/hgTables
#as a Gzipped .csv file, by setting Position to chr1:1-248,956,422, then chr2:1-242,193,529, etc.
#Both 'GH Reg Elems (DE)' and 'GH Interactions (DE)' tables were used.
#Saved to folder [workdir]/GH
mkdir(paste0(workdir,"/GH"))

#Import GH Reg Elems (DE) data for each chromosome and join:
library(data.table)
c1 <- fread(paste0(workdir,"/GH/c1.csv.gz"))
c2 <- fread(paste0(workdir,"/GH/c2.csv.gz"))
c3 <- fread(paste0(workdir,"/GH/c3.csv.gz"))
c4 <- fread(paste0(workdir,"/GH/c4.csv.gz"))
c5 <- fread(paste0(workdir,"/GH/c5.csv.gz"))
c6 <- fread(paste0(workdir,"/GH/c6.csv.gz"))
c7 <- fread(paste0(workdir,"/GH/c7.csv.gz"))
c8 <- fread(paste0(workdir,"/GH/c8.csv.gz"))
c9 <- fread(paste0(workdir,"/GH/c9.csv.gz"))
c10 <- fread(paste0(workdir,"/GH/c10.csv.gz"))
c11 <- fread(paste0(workdir,"/GH/c11.csv.gz"))
c12 <- fread(paste0(workdir,"/GH/c12.csv.gz"))
c13 <- fread(paste0(workdir,"/GH/c13.csv.gz"))
c14 <- fread(paste0(workdir,"/GH/c14.csv.gz"))
c15 <- fread(paste0(workdir,"/GH/c15.csv.gz"))
c16 <- fread(paste0(workdir,"/GH/c16.csv.gz"))
c17 <- fread(paste0(workdir,"/GH/c17.csv.gz"))
c18 <- fread(paste0(workdir,"/GH/c18.csv.gz"))
c19 <- fread(paste0(workdir,"/GH/c19.csv.gz"))
c20 <- fread(paste0(workdir,"/GH/c20.csv.gz"))
c21 <- fread(paste0(workdir,"/GH/c21.csv.gz"))
c22 <- fread(paste0(workdir,"/GH/c22.csv.gz"))
cX <- fread(paste0(workdir,"/GH/cX.csv.gz"))
cY <- fread(paste0(workdir,"/GH/cY.csv.gz"))

c_list <- list(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13,
               c14, c15, c16, c17, c18, c19, c20, c21, c22, cX, cY)
e <- rbindlist(c_list)
fwrite(e, paste0(workdir,'/GH/genehancers.csv.gz'))

#Import GH Interactions (DE) data for each chromosome and join:
ci1 <- fread(paste0(workdir,"/GH/ci1.csv.gz"))
ci2 <- fread(paste0(workdir,"/GH/ci2.csv.gz"))
ci3 <- fread(paste0(workdir,"/GH/ci3.csv.gz"))
ci4 <- fread(paste0(workdir,"/GH/ci4.csv.gz"))
ci5 <- fread(paste0(workdir,"/GH/ci5.csv.gz"))
ci6 <- fread(paste0(workdir,"/GH/ci6.csv.gz"))
ci7 <- fread(paste0(workdir,"/GH/ci7.csv.gz"))
ci8 <- fread(paste0(workdir,"/GH/ci8.csv.gz"))
ci9 <- fread(paste0(workdir,"/GH/ci9.csv.gz"))
ci10 <- fread(paste0(workdir,"/GH/ci10.csv.gz"))
ci11 <- fread(paste0(workdir,"/GH/ci11.csv.gz"))
ci12 <- fread(paste0(workdir,"/GH/ci12.csv.gz"))
ci13 <- fread(paste0(workdir,"/GH/ci13.csv.gz"))
ci14 <- fread(paste0(workdir,"/GH/ci14.csv.gz"))
ci15 <- fread(paste0(workdir,"/GH/ci15.csv.gz"))
ci16 <- fread(paste0(workdir,"/GH/ci16.csv.gz"))
ci17 <- fread(paste0(workdir,"/GH/ci17.csv.gz"))
ci18 <- fread(paste0(workdir,"/GH/ci18.csv.gz"))
ci19 <- fread(paste0(workdir,"/GH/ci19.csv.gz"))
ci20 <- fread(paste0(workdir,"/GH/ci20.csv.gz"))
ci21 <- fread(paste0(workdir,"/GH/ci21.csv.gz"))
ci22 <- fread(paste0(workdir,"/GH/ci22.csv.gz"))
ciX <- fread(paste0(workdir,"/GH/ciX.csv.gz"))
ciY <- fread(paste0(workdir,"/GH/ciY.csv.gz"))

ci_list <- list(ci1, ci2, ci3, ci4, ci5, ci6, ci7, ci8, ci9, ci10, ci11, ci12, ci13,
               ci14, ci15, ci16, ci17, ci18, ci19, ci20, ci21, ci22, ciX, ciY)
ei <- rbindlist(ci_list)
fwrite(ei, paste0(workdir,'/GH/genehancers_interactions.csv.gz'))

#=========================================
#11B. Get info for each GeneHancer element
#=========================================
#Including:
#- Element ID
#- Element position (chromosome, start, and end)
#- Element type (promoter, enhancer, or both)
#- Gene(s) regulated by the element
#- Methods of associating each gene with the element

library(dplyr)

ei$geneAssociationMethods <- gsub(',','&',ei$geneAssociationMethods) #replace commas to avoid csv errors

#Reformat GH Interactions (DE) data so there is one row per regulatory element (so multiple genes per row,
#if the the element regulates multiple genes), rather than one row per gene/element interaction:
ei2 <- ei %>% 
  group_by(geneHancerIdentifier) %>%
  summarise(geneName = paste(unique(geneName), collapse=';'),
            geneAssociationMethods = paste(unique(geneAssociationMethods), collapse=';'))

#Join relevant info from GH Interactions (DE) and GH Reg Elems (DE)
colnames(ei2)[1] <- 'name'
e2 <- e[,c(1:4,11)]
ea <- left_join(e2,ei2,by='name')
colnames(ea) <- c('CHR','GeneHancer_Start','GeneHancer_End','GeneHancer_Name','GeneHancer_Feature_Type','GeneHancer_Associated_Gene','GeneHancer_Association_Methods')

#===============================================
#11C. Temporarily join GH annotation to manifest
#===============================================
#Manifest already imported and set to data.table in step 8
#Extract CpG ID, position and chromosome
m3 <- manifest[,.(IlmnID,CHR,MAPINFO)]
#Set GH data to data.table
setDT(ea)
#Initialise m3 with column names from ea (empty columns)
colz <- colnames(ea)[colnames(ea) != 'CHR']
m3[, (colz):=character()]

save(ea, m3, paste0(workdir,'/GH_annotation_files.RData'))

#Split by chromosome again (saves computation time)
ea1 <- ea[CHR == 'chr1',]
ea2 <- ea[CHR == 'chr2',]
ea3 <- ea[CHR == 'chr3',]
ea4 <- ea[CHR == 'chr4',]
ea5 <- ea[CHR == 'chr5',]
ea6 <- ea[CHR == 'chr6',]
ea7 <- ea[CHR == 'chr7',]
ea8 <- ea[CHR == 'chr8',]
ea9 <- ea[CHR == 'chr9',]
ea10 <- ea[CHR == 'chr10',]
ea11 <- ea[CHR == 'chr11',]
ea12 <- ea[CHR == 'chr12',]
ea13 <- ea[CHR == 'chr13',]
ea14 <- ea[CHR == 'chr14',]
ea15 <- ea[CHR == 'chr15',]
ea16 <- ea[CHR == 'chr16',]
ea17 <- ea[CHR == 'chr17',]
ea18 <- ea[CHR == 'chr18',]
ea19 <- ea[CHR == 'chr19',]
ea20 <- ea[CHR == 'chr20',]
ea21 <- ea[CHR == 'chr21',]
ea22 <- ea[CHR == 'chr22',]
eaX <- ea[CHR == 'chrX',]
eaY <- ea[CHR == 'chrY',]
eaM <- ea[CHR == 'chrM',]

mf1 <- m3[CHR == 'chr1',]
mf2 <- m3[CHR == 'chr2',]
mf3 <- m3[CHR == 'chr3',]
mf4 <- m3[CHR == 'chr4',]
mf5 <- m3[CHR == 'chr5',]
mf6 <- m3[CHR == 'chr6',]
mf7 <- m3[CHR == 'chr7',]
mf8 <- m3[CHR == 'chr8',]
mf9 <- m3[CHR == 'chr9',]
mf10 <- m3[CHR == 'chr10',]
mf11 <- m3[CHR == 'chr11',]
mf12 <- m3[CHR == 'chr12',]
mf13 <- m3[CHR == 'chr13',]
mf14 <- m3[CHR == 'chr14',]
mf15 <- m3[CHR == 'chr15',]
mf16 <- m3[CHR == 'chr16',]
mf17 <- m3[CHR == 'chr17',]
mf18 <- m3[CHR == 'chr18',]
mf19 <- m3[CHR == 'chr19',]
mf20 <- m3[CHR == 'chr20',]
mf21 <- m3[CHR == 'chr21',]
mf22 <- m3[CHR == 'chr22',]
mfX <- m3[CHR == 'chrX',]
mfY <- m3[CHR == 'chrY',]
mfM <- m3[CHR == 'chrM',]

#==== Label which GH elements each CpG site is in based off start/end of element and MAPINFO of CpG =====

mf1_2 <- mf1 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea1[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf2_2 <- mf2 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea2[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf3_2 <- mf3 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea3[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf4_2 <- mf4 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea4[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf5_2 <- mf5 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea5[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf6_2 <- mf6 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea6[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf7_2 <- mf7 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea7[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf8_2 <- mf8 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea8[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf9_2 <- mf9 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea9[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf10_2 <- mf10 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea10[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf11_2 <- mf11 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea11[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf12_2 <- mf12 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea12[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf13_2 <- mf13 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea13[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf14_2 <- mf14 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea14[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf15_2 <- mf15 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea15[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf16_2 <- mf16 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea16[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,
     mf16_2, file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf17_2 <- mf17 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea17[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,
     mf16_2, mf17_2, file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf18_2 <- mf18 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea18[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,
     mf16_2, mf17_2, mf18_2, file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf19_2 <- mf19 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea19[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,
     mf16_2, mf17_2, mf18_2, mf19_2, file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf20_2 <- mf20 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea20[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,
     mf16_2, mf17_2, mf18_2, mf19_2, mf20_2, file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf21_2 <- mf21 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea21[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,
     mf16_2, mf17_2, mf18_2, mf19_2, mf20_2, mf21_2, file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mf22_2 <- mf22 %>%
  group_by(MAPINFO) %>%
  do({
    summarise(ea22[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,
     mf16_2, mf17_2, mf18_2, mf19_2, mf20_2, mf21_2,
     mf22_2, file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mfX_2 <- mfX %>%
  group_by(MAPINFO) %>%
  do({
    summarise(eaX[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,
     mf16_2, mf17_2, mf18_2, mf19_2, mf20_2, mf21_2,
     mf22_2, mfX_2, file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mfY_2 <- mfY %>%
  group_by(MAPINFO) %>%
  do({
    summarise(eaY[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,
     mf16_2, mf17_2, mf18_2, mf19_2, mf20_2, mf21_2,
     mf22_2, mfX_2, mfY_2, file=paste0(workdir,'/GH_manifest_by_chr.RData'))
mfM_2 <- mfM %>%
  group_by(MAPINFO) %>%
  do({
    summarise(eaM[GeneHancer_Start <= .$MAPINFO & GeneHancer_End >= .$MAPINFO,],
              across(everything(), ~ paste(., collapse = ";")))
  })
save(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,
     mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,
     mf16_2, mf17_2, mf18_2, mf19_2, mf20_2, mf21_2,
     mf22_2, mfX_2, mfY_2, mfM_2, file=paste0(workdir,'/GH_manifest_by_chr.RData'))

#======== recombine across chromosomes =================

mf_list <- list(mf1_2,mf2_2,mf3_2,mf4_2,mf5_2,mf6_2,mf7_2,mf8_2,mf9_2,mf10_2,mf11_2,mf12_2,mf13_2,mf14_2,mf15_2,mf16_2,mf17_2,mf18_2,mf19_2,mf20_2,mf21_2,mf22_2,mfX_2,mfY_2,mfM_2)
chrs <- c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12',
          'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX',
          'chrY','chrM')
for (i in 1:25) {
  mf_list[[i]]$CHR <- chrs[i]
}
m4 <- bind_rows(mf_list$mf1_2,mf_list$mf2_2,mf_list$mf3_2,mf_list$mf4_2,mf_list$mf5_2,mf_list$mf6_2,mf_list$mf7_2,mf_list$mf8_2,mf_list$mf9_2,mf_list$mf10_2,mf_list$mf11_2,mf_list$mf12_2,mf_list$mf13_2,mf_list$mf14_2,mf_list$mf15_2,mf_list$mf16_2,mf_list$mf17_2,mf_list$mf18_2,mf_list$mf19_2,mf_list$mf20_2,mf_list$mf21_2,mf_list$mf22_2,mf_list$mfX_2,mf_list$mfY_2,mf_list$mfM_2)
m4 <- mf_list[[1]]
for (i in 2:25) {
  m4 <- bind_rows(m4,mf_list[[i]])
}

save(m4,file=paste0(workdir,'/2_GH_manifest_by_chr_combined.RData'))

#========================================================================================
#11D. Add binary TRUE/FALSE indicator of whether each CpG site is in a GeneHancer element
#========================================================================================
m4$In_GeneHancer <- ifelse(m4$GeneHancer_Feature_Type == "", FALSE, TRUE)

save(m4,file=paste0(workdir,'/2_GH_manifest_by_chr_combined.RData'))

#=====================================================================================================
#12. Append GeneHancer annotation to reannotated manifest (currently containing GENCODEv49 annotation)
#=====================================================================================================
load(paste0(workdir,'/1_GENCODE_manifest_by_chr_combined.RData'))
load(paste0(workdir,'/2_GENCODE_manifest_by_chr_combined.RData'))

m2 <- as.data.frame(m2)
m4$CM <- paste0(m4$CHR,'_',m4$MAPINFO)

#control probes and those mapped to "chromosome 0" should have no annotation
controls <- m2[m2$CM %in% c('chr0_0','_NA'),]
controls[,c("GENCODEv49_Gene_Name","GENCODEv49_Gene_ID",
                  "GENCODEv49_Gene_Type","GENCODEv49_Transcript_Name",
                  "GENCODEv49_Transcript_ID","GENCODEv49_Feature_Type",
                  "GENCODEv49_Feature_Exon_Number","GENCODEv49_Feature_Type_Specific",
                  "GENCODEv49_In_Gene_Body","GENCODEv49_Gene_Start",
                  "GENCODEv49_Gene_End","GENCODEv49_Gene_Strand",
                  "GENCODEv49_Transcript_Start","GENCODEv49_Transcript_End",
                  "GENCODEv49_Feature_Start","GENCODEv49_Feature_End",
                  "In_GeneHancer","Horvath_Multiple_Clock_Site","MethylDetectR_Clock_Site",
                  "GeneHancer_Start","GeneHancer_End","GeneHancer_Name","GeneHancer_Feature_Type",
                  "GeneHancer_Associated_Gene","GeneHancer_Association_Methods")] <- NA
controls <- controls[,colnames(controls) != 'CM']

#left join GENCODE and GeneHancer annotations, excluding controls/chr0 (which will interfere with left join)
m5 <- dplyr::left_join(m2[!(m2$CM %in% c('chr0_0','_NA')),],
                       m4[!(m4$CM %in% c('chr0_0','_NA')),c('CM','In_GeneHancer',
                                 'GeneHancer_Start','GeneHancer_End','GeneHancer_Name',
                                 'GeneHancer_Feature_Type','GeneHancer_Associated_Gene',
                                 'GeneHancer_Association_Methods')], by='CM', keep=FALSE)

#=======================================================================
#11F. Remove all GeneHancer data from manifest (as can't share/reupload)
#=======================================================================
m6 <- m5[,!(colnames(m5) %in% c('GeneHancer_Start','GeneHancer_End',
                                'GeneHancer_Name','GeneHancer_Feature_Type',
                                'GeneHancer_Associated_Gene',
                                'GeneHancer_Association_Methods'))]

save(m6,file=paste0(workdir,'/3_GENCODE_public_manifest.RData'))

#=========================================================
#12. Add information from Horvath and MethylDetectR clocks
#=========================================================
#load in csvs of CpG site names from Horvath and MethylDetectR
h <- read.csv(paste0(workdir,'/clocks/datMiniAnnotation4_fixed.csv')) #Horvath clock sites downloaded from https://dnamage.clockfoundation.org/importanthints-page 
mdr <- read.csv(paste0(workdir,'/clocks/Truncate_to_these_CpGs.csv')) #MethylDetectR clock sites downloaded from https://zenodo.org/records/7154750 

#add to manifest based on EPICv1_Loci column
m6$Horvath_Multiple_Clock_Site <- ifelse(m6$EPICv1_Loci %in% h$Name, TRUE, FALSE)
m6$MethylDetectR_Clock_Site <- ifelse(m6$EPICv1_Loci %in% mdr$CpGs, TRUE, FALSE)

save(m6,file=paste0(workdir,'/4_GENCODE_GH_clocks_public_manifest.RData'))

#==============================
#13. Add back in control probes
#==============================
m7 <- m6[,colnames(m6) != 'CM']

unique((colnames(m7) %in% colnames(controls))) #TRUE]
controls7 <- controls[,colnames(controls) %in% colnames(m7)]
identical(colnames(m7), colnames(controls7)) #TRUE

m8 <- rbind(m7, controls7)

#reorder into same order as original manifest
probe_order <- data.frame(IlmnID=manifest$IlmnID)
m9 <- left_join(probe_order, m8)
identical(m9$IlmnID, manifest$IlmnID) #TRUE

#=====================================================================
#14. Remove redundant Illumina columns (e.g. from GENCODEv41 and UCSC)
#=====================================================================
m10 <- m9[,!(colnames(m9) %in% c('UCSC_RefGene_Group','UCSC_RefGene_Name','UCSC_RefGene_Accession',
                'GencodeV41_Group','GencodeV41_Name','GencodeV41_Accession','CM'))]

#===============================
#15. Export for upload to Zenodo
#===============================
system(paste0('mkdir ',workdir,'/FINAL'))

save(m9, file=paste0(workdir,'/FINAL/EPICv2_reannotated_manifest_forcomparison_v3.0.RData'))
save(m10, file=paste0(workdir,'/FINAL/EPICv2_reannotated_manifest_v3.0.RData'))
fwrite(m9, paste0(workdir,'/FINAL/EPICv2_reannotated_manifest_forcomparison_v3.0.csv.gz'), row.names=FALSE)
fwrite(m10, paste0(workdir,'/FINAL/EPICv2_reannotated_manifest_v3.0.csv.gz'), row.names=FALSE)
