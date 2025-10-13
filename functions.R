require(tidyr)
require(dplyr)

expand_annotation <- function(manifest, by='gene', verbose=TRUE) {
  if (verbose) {
    if (by=='gene') {
      print('Returning CpGs with annotated GENCODEv47 genes.')
      print('Important! This function assumes your CpG IDs/IlmnIDs are in the first column of input file, and the input file contains columns named GENCODEv47_Gene_ID and/or GENCODEv47_Gene_Name.')
    }
    if (by=='transcript') {
      print('Returning CpGs with annotated GENCODEv47 transcripts.')
      print('Important! This function assumes your CpG IDs/IlmnIDs are in the first column of input file, and the input file contains columns named GENCODEv47_Transcript_ID and/or GENCODEv47_Transcript_Name.')
    }
    if (by=='db') {
      print('Returning CpGs with annotated regulatory elements (based on distance from genes).')
      print('Important! This function assumes your CpG IDs/IlmnIDs are in the first column of input file, and the input file contains columns named DB_Element_Gene_ID and/or DB_Element_Gene_Name')
    }
    print('Set verbose = FALSE to silence this message.')
  }

  if (by=='gene') {
    g47_cols <- colnames(manifest)[grepl('GENCODEv47',colnames(manifest))]
    gname <- 'GENCODEv47_Gene_Name' %in% g47_cols
    gid  <- 'GENCODEv47_Gene_ID' %in% g47_cols

  
    if (!gname & !gid) {
      stop('Columns GENCODEv47_Gene_Name and GENCODEv47_Gene_ID not found in input file.')
    }
  
    temp <- tidyr::separate_longer_delim(manifest,
                                       cols=all_of(g47_cols),
                                       delim=';')
  
    if (gid) {
      temp$temp <- paste0(temp[,1], temp$GENCODEv47_Gene_ID)
    } else {
      temp$temp <- paste0(temp[,1], temp$GENCODEv47_Gene_Name)
    }
  
    temp <- temp[!duplicated(temp$temp),]
  
    if (gname) {
      mgenes <- temp[,colnames(temp) %in% c(colnames(temp)[1], 'GENCODEv47_Gene_Name', 'temp')]
    } else {
      mgenes <- temp[,colnames(temp) %in% c(colnames(temp)[1], 'GENCODEv47_Gene_ID', 'temp')]
    }
    colnames(mgenes)[2] <- 'Gene'
    mgenes <- left_join(mgenes, temp[,-1], by='temp')
    mgenes <- mgenes[,which(colnames(mgenes) != 'temp')]
    return(mgenes)
  }
  if (by=='transcript') {
    g47_cols <- colnames(manifest)[grepl('GENCODEv47',colnames(manifest))]
    tname <- 'GENCODEv47_Transcript_Name' %in% g47_cols
    tid  <- 'GENCODEv47_Transcript_ID' %in% g47_cols

  
    if (!tname & !tid) {
      stop('Columns GENCODEv47_Transcript_Name and GENCODEv47_Transcript_ID not found in input file.')
    }
  
    temp <- tidyr::separate_longer_delim(manifest,
                                       cols=all_of(g47_cols),
                                       delim=';')
  
    if (tid) {
      temp$temp <- paste0(temp[,1], temp$GENCODEv47_Transcript_ID)
    } else {
      temp$temp <- paste0(temp[,1], temp$GENCODEv47_Transcript_Name)
    }
  
    temp <- temp[!duplicated(temp$temp),]
  
    if (tname) {
      mtrans <- temp[,colnames(temp) %in% c(colnames(temp)[1], 'GENCODEv47_Transcript_Name', 'temp')]
    } else {
      mtrans <- temp[,colnames(temp) %in% c(colnames(temp)[1], 'GENCODEv47_Transcript_ID', 'temp')]
    }
    colnames(mtrans)[2] <- 'Transcript'
    mtrans <- left_join(mtrans, temp[,-1], by='temp')
    mtrans <- mtrans[,which(colnames(mtrans) != 'temp')]
    return(mtrans)
  }
  if (by=='db') {
    db_cols <- colnames(manifest)[grepl('DB_Element',colnames(manifest))]
    gname <- 'DB_Element_Gene_Name' %in% db_cols
    gid  <- 'DB_Element_Gene_ID' %in% db_cols
    etype  <- 'DB_Element_Type' %in% db_cols

  
    if (!gname & !gid) {
      stop('Columns DB_Element_Gene_Name and DB_Element_Gene_ID not found in input file.')
    }
  
    temp <- tidyr::separate_longer_delim(manifest,
                                       cols=all_of(db_cols),
                                       delim=';')
  
    if (gid) {
      temp$temp <- paste0(temp[,1], temp$DB_Element_Gene_ID)
    } else {
      temp$temp <- paste0(temp[,1], temp$DB_Element_Gene_Name)
    }
  
    temp <- temp[!duplicated(temp$temp),]
  
    if (gname) {
      mgenes <- temp[,colnames(temp) %in% c(colnames(temp)[1], 'DB_Element_Gene_Name', 'temp')]
    } else {
      mgenes <- temp[,colnames(temp) %in% c(colnames(temp)[1], 'DB_Element_Gene_ID', 'temp')]
    }
    if (etype) {
      mgenes$Element <- rep('',nrow(mgenes))
      mgenes$Element[grepl('Promoter',temp$DB_Element_Type,ignore.case = TRUE)] <- 'Promoter'
      mgenes$Element[grepl('Enhancer',temp$DB_Element_Type,ignore.case = TRUE)] <- 'Enhancer'
    }
    colnames(mgenes)[2] <- 'Gene'
    mgenes <- left_join(mgenes, temp[,-1], by='temp')
    mgenes <- mgenes[,which(colnames(mgenes) != 'temp')]
    return(mgenes)
  }
  
}
