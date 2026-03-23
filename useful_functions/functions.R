require(tidyr)
require(dplyr)

expand_annotation <- function(manifest, by='gene', verbose=TRUE) {
  if (verbose) {
    if (by=='gene') {
      print('Returning CpGs with annotated genes from GENCODEv47 (if using reannotated manifest version 2.0 or older) or GENCODEv49 (version 3.0+).')
      print('Important! This function assumes your CpG IDs/IlmnIDs are in the first column of input file, and the input file contains at least one of the following columns: GENCODEv47_Gene_ID, GENCODEv47_Gene_Name, GENCODEv49_Gene_ID, GENCODEv49_Gene_Name')
    }
    if (by=='transcript') {
      print('Returning CpGs with annotated transcripts from GENCODEv47 (if using reannotated manifest version 2.0 or older) or GENCODEv49 (version 3.0+).')
      print('Important! This function assumes your CpG IDs/IlmnIDs are in the first column of input file, and the input file contains at least one of the following columns: GENCODEv47_Transcript_ID, GENCODEv47_Transcript_Name, GENCODEv49_Transcript_ID, GENCODEv49_Transcript_Name')
    }
    if (by=='gh') {
      print('Returning CpGs annotated to regulatory elements (according to GeneHancer database).')
      print('Important! This function assumes your CpG IDs/IlmnIDs are in the first column of input file, and the input file contains columns named In_GeneHancer and/or GeneHancer_Associated_Gene.')
      print('If you are using the publicly available version of the manifest, you will only have column In_GeneHancer, due to GeneHancer T&Cs.')
    }
    print('Set verbose = FALSE to silence this message.')
  }

  if (by=='gene') {
    g_cols <- colnames(manifest)[grepl('GENCODEv47|GENCODEv49',colnames(manifest))]
    g47name <- 'GENCODEv47_Gene_Name' %in% g_cols
    g47id  <- 'GENCODEv47_Gene_ID' %in% g_cols
    g49name <- 'GENCODEv49_Gene_Name' %in% g_cols
    g49id  <- 'GENCODEv49_Gene_ID' %in% g_cols
    gname <- g47name | g49name
    gid <- g47id | g49id

  
    if (!gname & !gid) {
      stop('At least one of the columns GENCODEv47_Gene_Name / GENCODEv47_Gene_ID / GENCODEv49_Gene_Name / GENCODEv49_Gene_ID is required in input file.')
    }

    if (grepl('v47', g_cols) & grepl('v49', g_cols)) {
      stop("This error means you have both the GENCODEv47 and GENCODEv49 annotations in your input file. Did you join the older and newer versions of the reannotated EPICv2 manifest together for some reason? I'm not judging, just confused.")
    }
  
    temp <- tidyr::separate_longer_delim(manifest,
                                       cols=all_of(g_cols),
                                       delim=';')
  
    if (g47id) {
      temp$temp <- paste0(temp[,1], temp$GENCODEv47_Gene_ID)
      gcol <- 'GENCODEv47_Gene_ID'
    } else if (g47name) {
      temp$temp <- paste0(temp[,1], temp$GENCODEv47_Gene_Name)
      gcol <- 'GENCODEv47_Gene_Name'
    } else if (g49id) {
      temp$temp <- paste0(temp[,1], temp$GENCODEv49_Gene_ID)
      gcol <- 'GENCODEv49_Gene_ID'
    } else if (g49name) {
      temp$temp <- paste0(temp[,1], temp$GENCODEv49_Gene_Name)
      gcol <- 'GENCODEv49_Gene_Name'
    }
  
    temp <- temp[!duplicated(temp$temp),]
  
    mgenes <- temp[,colnames(temp) %in% c(colnames(temp)[1], gcol, 'temp')]
    colnames(mgenes)[2] <- 'Gene'
    
    mgenes <- left_join(mgenes, temp[,-1], by='temp')
    mgenes <- mgenes[,which(colnames(mgenes) != 'temp')]
    return(mgenes)
  }
  if (by=='transcript') {
    g_cols <- colnames(manifest)[grepl('GENCODEv47|GENCODEv49',colnames(manifest))]
    g47name <- 'GENCODEv47_Transcript_Name' %in% g_cols
    g47id  <- 'GENCODEv47_Transcript_ID' %in% g_cols
    g49name <- 'GENCODEv49_Transcript_Name' %in% g_cols
    g49id  <- 'GENCODEv49_Transcript_ID' %in% g_cols
    gname <- g47name | g49name
    gid <- g47id | g49id

  
    if (!gname & !gid) {
      stop('At least one of the columns GENCODEv47_Transcript_Name / GENCODEv47_Transcript_ID / GENCODEv49_Transcript_Name / GENCODEv49_Transcript_ID is required in input file.')
    }

    if (grepl('v47', g_cols) & grepl('v49', g_cols)) {
      stop("This error means you have both the GENCODEv47 and GENCODEv49 annotations in your input file. Did you join the older and newer versions of the reannotated EPICv2 manifest together for some reason? (Not judging, just confused.)")
    }
  
    temp <- tidyr::separate_longer_delim(manifest,
                                       cols=all_of(g_cols),
                                       delim=';')
  
    if (g47id) {
      temp$temp <- paste0(temp[,1], temp$GENCODEv47_Transcript_ID)
      gcol <- 'GENCODEv47_Transcript_ID'
    } else if (g47name) {
      temp$temp <- paste0(temp[,1], temp$GENCODEv47_Transcript_Name)
      gcol <- 'GENCODEv47_Transcript_Name'
    } else if (g49id) {
      temp$temp <- paste0(temp[,1], temp$GENCODEv49_Transcript_ID)
      gcol <- 'GENCODEv49_Transcript_ID'
    } else if (g49name) {
      temp$temp <- paste0(temp[,1], temp$GENCODEv49_Transcript_Name)
      gcol <- 'GENCODEv49_Transcript_Name'
    }
  
    temp <- temp[!duplicated(temp$temp),]
  
    mtrans <- temp[,colnames(temp) %in% c(colnames(temp)[1], gcol, 'temp')]
    colnames(mtrans)[2] <- 'Transcript'
    
    mtrans <- left_join(mtrans, temp[,-1], by='temp')
    mtrans <- mtrans[,which(colnames(mtrans) != 'temp')]
    return(mtrans)
  }
  
  if (by=='gh') {
    gh_cols <- colnames(manifest)[grepl('GeneHancer',colnames(manifest))]
    ig <- 'In_GeneHancer' %in% gh_cols
    gname <- 'GeneHancer_Associated_Gene' %in% gh_cols
    ghname <- 'GeneHancer_Name' %in% gh_cols
    etype  <- 'GeneHancer_Feature_Type' %in% gh_cols
    print('debug 1')
    
    if(!ig & !gname) {
      stop('Columns In_GeneHancer and GeneHancer_Associated_Gene not found in input file.')
    }
      
    if(ig & !gname & !ghname & !etype) {
      print('debug 2')
      colorder <- c(1,which(colnames(manifest)=='In_GeneHancer'),(2:ncol(manifest))[-(which(colnames(manifest)=='In_GeneHancer')-1)])
      return(manifest[,..colorder])
    } else {
      print('debug 3')
      temp <- tidyr::separate_longer_delim(manifest,
                                       cols='GeneHancer_Associated_Gene',
                                       delim=';')
      print('debug 4')
      
      if (gname & ghname) {
        temp$temp <- paste0(temp[,1], temp$GeneHancer_Associated_Gene, temp$GeneHancer_Name)
        temp <- temp[!duplicated(temp$temp),]
        if (etype) {
          mgh <- data.frame(CpG = temp[,1], Gene = temp$GeneHancer_Associated_Gene, Element = temp$GeneHancer_Feature_Type,
                           GeneHancer_Element_Name = temp$GeneHancer_Name, temp = temp$temp)
          colnames(mgh)[1] <- colnames(temp)[1]
        } else {
          mgh <- data.frame(CpG = temp[,1], Gene = temp$GeneHancer_Associated_Gene,
                           GeneHancer_Element_Name = temp$GeneHancer_Name, temp = temp$temp)
          colnames(mgh)[1] <- colnames(temp)[1]
        }
      } else if (gname & !ghname) {
        temp$temp <- paste0(temp[,1], temp$GeneHancer_Associated_Gene)
        temp <- temp[!duplicated(temp$temp),]
        if (etype) {
          mgh <- data.frame(CpG = temp[,1], Gene = temp$GeneHancer_Associated_Gene, Element = temp$GeneHancer_Feature_Type, temp = temp$temp)
          colnames(mgh)[1] <- colnames(temp)[1]
        } else {
          mgh <- data.frame(CpG = temp[,1], Gene = temp$GeneHancer_Associated_Gene, temp = temp$temp)
          colnames(mgh)[1] <- colnames(temp)[1]
        }
      } else if (ghname & !gname) {
        temp$temp <- paste0(temp[,1], temp$GeneHancer_Name)
        temp <- temp[!duplicated(temp$temp),]
        if (etype) {
          mgh <- data.frame(CpG = temp[,1], Element = temp$GeneHancer_Feature_Type,
                           GeneHancer_Element_Name = temp$GeneHancer_Name, temp = temp$temp)
          colnames(mgh)[1] <- colnames(temp)[1]
        } else {
          mgh <- data.frame(CpG = temp[,1], GeneHancer_Element_Name = temp$GeneHancer_Name, temp = temp$temp)
          colnames(mgh)[1] <- colnames(temp)[1]
        }
      }
      mgh <- left_join(mgh, temp[,-1], by='temp')
      mgh <- mgh[,which(colnames(mgh) != 'temp')]
      return(mgh)
    }
  }
  
}

filter_to_genebody <- function(expanded_annotation, sig_cpgs = NULL) {
  if (sum(grepl('GENCODEv.._Feature_Type',colnames(expanded_annotation)))) {
    stop('Column named GENCODEv47_Feature_Type or GENCODEv49_Feature_Type missing from input file.')
  }
  gcol <- colnames(expanded_annotation)[grepl('GENCODEv.._Feature_Type',colnames(expanded_annotation))]
  if (sum(grepl(';', expanded_annotation[,gcol])) > 0) {
    stop("Did you apply the function expand_annotation() with by='gene' to your input file before running this function?")
  }
  if (!is.null(sig_cpgs)) {
    temp <- expanded_annotation[expanded_annotation[,1] %in% sig_cpgs,]
    temp <- temp[grepl('intron|exon|UTR',temp[,gcol]),]
  } else {
    temp <- expanded_annotation[grepl('intron|exon|UTR',expanded_annotation[,gcol]),]
  }
  return(temp)
}

filter_to_TSS1500 <- function(expanded_annotation, sig_cpgs = NULL) {
  if (!('GENCODEv47_Feature_Type' %in% colnames(expanded_annotation))) {
    stop('Column named GENCODEv47_Feature_Type missing from input file.')
  }
  if (sum(grepl(';', expanded_annotation$GENCODEv47_Feature_Type)) > 0) {
    stop("Did you apply the function expand_annotation() with by='gene' to your input file before running this function?")
  }
  if (!is.null(sig_cpgs)) {
    temp <- expanded_annotation[expanded_annotation[,1] %in% sig_cpgs,]
    temp <- temp[grepl('TSS1500',temp$GENCODEv47_Feature_Type),]
  } else {
    temp <- expanded_annotation[grepl('TSS1500',expanded_annotation$GENCODEv47_Feature_Type),]
  }
  return(temp)
}

filter_to_TSS200 <- function(expanded_annotation, sig_cpgs = NULL) {
  if (!('GENCODEv47_Feature_Type' %in% colnames(expanded_annotation))) {
    stop('Column named GENCODEv47_Feature_Type missing from input file.')
  }
  if (sum(grepl(';', expanded_annotation$GENCODEv47_Feature_Type)) > 0) {
    stop("Did you apply the function expand_annotation() with by='gene' to your input file before running this function?")
  }
  if (!is.null(sig_cpgs)) {
    temp <- expanded_annotation[expanded_annotation[,1] %in% sig_cpgs,]
    temp <- temp[grepl('TSS200',temp$GENCODEv47_Feature_Type),]
  } else {
    temp <- expanded_annotation[grepl('TSS200',expanded_annotation$GENCODEv47_Feature_Type),]
  }
  return(temp)
}

filter_to_promoter <- function(expanded_annotation, sig_cpgs = NULL) {
  if (!('DB_Element_Type' %in% colnames(expanded_annotation))) {
    stop('Column named DB_Element_Type missing from input file.')
  }
  if (sum(grepl(';', expanded_annotation$DB_Element_Type)) > 0) {
    stop("Did you apply the function expand_annotation() with by='db' to your input file before running this function?")
  }
  if (!is.null(sig_cpgs)) {
    temp <- expanded_annotation[expanded_annotation[,1] %in% sig_cpgs,]
    temp <- temp[grepl('Promoter',temp$DB_Element_Type,ignore.case=TRUE),]
  } else {
    temp <- expanded_annotation[grepl('Promoter',expanded_annotation$DB_Element_Type,ignore.case=TRUE),]
  }
  return(temp)
}

filter_to_enhancer <- function(expanded_annotation, sig_cpgs = NULL) {
  if (!('DB_Element_Type' %in% colnames(expanded_annotation))) {
    stop('Column named DB_Element_Type missing from input file.')
  }
  if (sum(grepl(';', expanded_annotation$DB_Element_Type)) > 0) {
    stop("Did you apply the function expand_annotation() with by='db' to your input file before running this function?")
  }
  if (!is.null(sig_cpgs)) {
    temp <- expanded_annotation[expanded_annotation[,1] %in% sig_cpgs,]
    temp <- temp[grepl('Enhancer',temp$DB_Element_Type,ignore.case=TRUE),]
  } else {
    temp <- expanded_annotation[grepl('Enhancer',expanded_annotation$DB_Element_Type,ignore.case=TRUE),]
  }
  return(temp)
}

filter_to_gh_promoter <- function(expanded_annotation, sig_cpgs = NULL) {
  if (!('GeneHancer_Feature_Type' %in% colnames(expanded_annotation))) {
    stop('Column named GeneHancer_Feature_Type missing from input file.')
  }
  if (sum(grepl(';', expanded_annotation$GeneHancer_Feature_Type)) > 0) {
    stop("Did you apply the function expand_annotation() with by='gh' to your input file before running this function?")
  }
  if (!is.null(sig_cpgs)) {
    temp <- expanded_annotation[expanded_annotation[,1] %in% sig_cpgs,]
    temp <- temp[grepl('Promoter',temp$GeneHancer_Feature_Type,ignore.case=TRUE),]
  } else {
    temp <- expanded_annotation[grepl('Promoter',expanded_annotation$GeneHancer_Feature_Type,ignore.case=TRUE),]
  }
  return(temp)
}

filter_to_gh_enhancer <- function(expanded_annotation, sig_cpgs = NULL) {
  if (!('GeneHancer_Feature_Type' %in% colnames(expanded_annotation))) {
    stop('Column named GeneHancer_Feature_Type missing from input file.')
  }
  if (sum(grepl(';', expanded_annotation$GeneHancer_Feature_Type)) > 0) {
    stop("Did you apply the function expand_annotation() with by='gh' to your input file before running this function?")
  }
  if (!is.null(sig_cpgs)) {
    temp <- expanded_annotation[expanded_annotation[,1] %in% sig_cpgs,]
    temp <- temp[grepl('Enhancer',temp$GeneHancer_Feature_Type,ignore.case=TRUE),]
  } else {
    temp <- expanded_annotation[grepl('Enhancer',expanded_annotation$GeneHancer_Feature_Type,ignore.case=TRUE),]
  }
  return(temp)
}

get_annotated_gene_list <- function(expanded_annotation, na.rm=FALSE) {
  if (na.rm) {
    return(unique(expanded_annotation$Gene[!is.na(expanded_annotation$Gene) & expanded_annotation$Gene != ""]))
  } else {
    return(unique(expanded_annotation$Gene))
  }
}

get_annotated_gene_table <- function(expanded_annotation, na.rm=FALSE) {
  if (na.rm) {
    return(table(expanded_annotation$Gene[!is.na(expanded_annotation$Gene) & expanded_annotation$Gene != ""]))
  } else {
    return(table(expanded_annotation$Gene))
  }
}

get_annotated_transcript_list <- function(expanded_annotation, na.rm=FALSE) {
  if (na.rm) {
    return(unique(expanded_annotation$Transcript[!is.na(expanded_annotation$Transcript) & expanded_annotation$Transcript != ""]))
  } else {
    return(unique(expanded_annotation$Transcript))
  }
}

get_annotated_transcript_table <- function(expanded_annotation, na.rm=FALSE) {
  if (na.rm) {
    return(table(expanded_annotation$Transcript[!is.na(expanded_annotation$Transcript) & expanded_annotation$Transcript != ""]))
  } else {
    return(table(expanded_annotation$Transcript))
  }
}
