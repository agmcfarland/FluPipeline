
tmpFile <- function(){ paste0('tmp.', paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = '')) }

samMD2length <- function(x){
  x <- sub('^MD\\:Z\\:', '', x)
  x <- gsub('\\^', '', x) # Deletions not supported.
  a <- as.integer(unlist(strsplit(x, '[A-Z]')))
  sum(a[!is.na(a)]) + length(unlist(str_match_all(x, '[A-Z]')))
}


str_diff <- function(x, y){
  o <- mapply(function(a, b){ 
    if(a != b) return(paste(a, '->', b))
  }, unlist(strsplit(x, '')), unlist(strsplit(y, ''))) 
  
  unlist(o[!sapply(o, is.null)])
  
}

shortRead2DNAstringSet <- function(x){
  r <- x@sread
  names(r) <- sub('\\s+.+$', '', as.character(x@id))
  r
}

prepareTrimmedReads <- function(R1, R2, qualCode = '5'){
  R1 <- trimTailw(R1, 2, qualCode, 5)
  R2 <- trimTailw(R2, 2, qualCode, 5)
  
  R1 <- shortRead2DNAstringSet(R1)
  R2 <- shortRead2DNAstringSet(R2)
  n  <- intersect(names(R1), names(R2))
  R1 <- R1[names(R1) %in% n]
  R2 <- R2[names(R2) %in% n]
  
  list(R1, R2)
}

createSoftwareVersionTable <- function(){
  if(! dir.exists(opt$samtoolsBin)) stop('Error - samtools bin does not exist')
  opt$samtoolsVersion <<- paste0(system(paste(file.path(opt$samtoolsBin, 'samtools'), '--version'), intern = TRUE)[1:2], collapse = ' ')
  opt$samtoolsVersion <<- sub('samtools\\s+', '', opt$samtoolsVersion)
  
  if(! dir.exists(opt$bcftoolsBin)) stop('Error - bcftools bin does not exist')
  opt$bcftoolsVersion <<- paste0(system(paste(file.path(opt$bcftoolsBin, 'bcftools'), '--version'), intern = TRUE)[1:2], collapse = ' ')
  opt$bcftoolsVersion <<- sub('bcftools\\s+', '', opt$bcftoolsVersion)
  
  if(! file.exists(opt$bwaPath)) stop('Error - bwa path does not exist')
  system(paste(opt$bwaPath, '2>', file.path(opt$workDir, 'bwa.version')))
  v <- readLines(file.path(opt$workDir, 'bwa.version'))
  opt$bwaVersion <<- sub('version:\\s+', '', tolower(v[grep('Version', v, ignore.case = TRUE)]))
  
  
  p <- c('#!/bin/bash', paste0('source ', opt$condaShellPath), 'conda activate pangolin', 'pangolin -v')
  writeLines(p, file.path(opt$workDir, 'pangolin.version.script'))
  system(paste0('chmod 755 ', file.path(opt$workDir, 'pangolin.version.script')))
  opt$pangolinVersion <<- sub('pangolin\\s+', '',  system(file.path(opt$workDir, 'pangolin.version.script'), intern = TRUE))
  
  p <- installed.packages()[names(sessionInfo()$otherPkgs), 'Version']
  bind_rows(tibble('Software/R package' = c('R', 'bwa', 'samtools', 'bcftools', 'pangolin'),
                   'Version' = c(paste0(R.Version()$major, '.', R.Version()$minor),
                                 opt$bwaVersion, opt$samtoolsVersion, opt$bcftoolsVersion, opt$pangolinVersion)),
            tibble('Software/R package' = names(p), 'Version' = p))
}


# megaHitContigs <- function(R1, R2, label = 'x', workDir = '.', megahit.path = 'megahit'){
#   contigs <- Reduce('append', lapply(c(1, 2, 3), function(n){
#     system(paste0(megahit.path, ' -t 15 -1 ', R1, ' -2 ', R2, '  -o ', workDir, 
#                   ' --min-count ', n, ' --k-list ', paste0(seq(from = 21, to = 127, by = 6), collapse = ',')))
#     
#     if(file.exists(file.path(workDir, 'final.contigs.fa'))){
#       contigs <- readFasta(file.path(workDir, 'final.contigs.fa'))
#     } else {
#       contigs <- ShortRead::ShortRead()
#     }
#     
#     unlink(workDir, recursive = TRUE)
#     contigs
#   }))
#   
#   contigs <- unique(contigs@sread)
#   contigs <- contigs[order(width(contigs), decreasing = TRUE)]
#   if(length(contigs) > 0) names(contigs) <- paste0('contig', 1:length(contigs))
#   contigs
# }


# Basic megaHit call.
megaHitContigs <- function(R1, R2, label = 'x', workDir = '.', megahit.path = 'megahit'){
  system(paste0(megahit.path, ' -t 15 -1 ', R1, ' -2 ', R2, '  -o ', workDir,
                ' --min-count 2 --k-list ', paste0(seq(from = 21, to = 127, by = 6), collapse = ',')))

  if(file.exists(file.path(workDir, 'final.contigs.fa'))){
    contigs <- readFasta(file.path(workDir, 'final.contigs.fa'))
    contigs <- unique(contigs@sread)
    contigs <- contigs[order(width(contigs), decreasing = TRUE)]
    if(length(contigs) > 0) names(contigs) <- paste0('contig', 1:length(contigs))
    return(contigs)
  } else {
    return(Biostrings::DNAStringSet())
  }
}



parsePileUpString <- function(x){
  # Remove read boundary markers
  x <- toupper(gsub('\\$|\\^\\S', '', x))
  
  # Remove leading strand indicator from indels. 
  x <- gsub('[\\.\\,]([\\-\\+])', '\\1', x, perl = TRUE)
  
  counts <- list('e' = 0, '?' = 0)
  
  indelStart <- FALSE
  indelLength <- 0
  indelCapture <- ''
  indelType <- NA
  
  for(i in unlist(strsplit(x, ''))){ 
    if(indelStart == TRUE){
      indelLength <- as.integer(i)
      indelStart <- FALSE
    } else if (indelLength > 1){
      indelCapture <- paste0(indelCapture, i)
      indelLength <- indelLength - 1
    }
    else if(indelLength == 1){
      indelCapture <- paste0(indelCapture, i)
      if(! paste0(indelType, indelCapture) %in% names(counts)) counts[[paste0(indelType, indelCapture)]] <- 0
      counts[[paste0(indelType, indelCapture)]] <- counts[[paste0(indelType, indelCapture)]] + 1
      indelCapture <- ''
      indelLength <- 0
    }
    else if (i == '.' || i == ','){
      counts$e <- counts$e + 1
    } else if (grepl('[ATCGN]', i)){
      if(! i %in% names(counts)) counts[[i]] <- 0
      counts[[i]] <- counts[[i]] + 1
    } else if (grepl('[\\+\\-]', i)){
      indelStart <- TRUE
      indelType <- ifelse(grepl('\\+', i), 'ins', 'del')
    } else {
      counts[['?']] <- counts[['?']] + 1
    }
  }
  
  unlist(counts)/sum(unlist(counts))
}


parsePileUpString2 <- function(x){
  # Remove read boundary markers and deletion markers.
  x <- toupper(gsub('\\*|\\$|\\^\\S', '', x))
  
  counts <- list('e' = 0, '?' = 0)
  
  indelStart <- FALSE
  indelLength <- ''
  indelCapture <- ''
  indelType <- NA
  indelCaptureCount <- 0
  
  for(i in unlist(strsplit(x, ''))){ 
    if(indelStart == TRUE & ! grepl('[ATCGN]', i)){  #  + or - was seen before and now we start to capture the length of the indel. 
        indelLength <- as.integer(paste0(indelLength, as.integer(i)))
    } else if (indelStart == TRUE & grepl('[ATCGN]', i)){  # The indel length has now been captured and now we record the first NT of the indel.
        indelStart <- FALSE
        indelCapture <- i
        indelCaptureCount <- 1
    } else if(indelCaptureCount < indelLength){   # Continue to capture the indel sequence.
        indelCapture <- paste0(indelCapture, i)
        indelCaptureCount <- indelCaptureCount + 1
    } else if(indelCaptureCount == indelLength){  # We've reached the end of the indel, record it.
        if(! paste0(indelType, indelCapture) %in% names(counts)) counts[[paste0(indelType, indelCapture)]] <- 0  # Set counter to zero if not seen.
        counts[[paste0(indelType, indelCapture)]] <- counts[[paste0(indelType, indelCapture)]] + 1
        indelLength  <- ''
        indelCapture <- ''
        indelCaptureCount <- 0
    } else if (i == '.' || i == ','){ # Encountered expected sequence of . or ,.
        counts$e <- counts$e + 1
    } else if (grepl('[ATCGN]', i)){  # Encountered vairant position.
        if(! i %in% names(counts)) counts[[i]] <- 0
        counts[[i]] <- counts[[i]] + 1
    } else if (grepl('[\\+\\-]', i)){ # Encounter a + or -, start recording an indel record on next loop starting with integer capture.
        indelStart <- TRUE
        indelType <- ifelse(grepl('\\+', i), 'ins', 'del')
    } else {                          # Fail safe, capture unknown events.
        counts[['?']] <- counts[['?']] + 1
    }
  }
  
  unlist(counts)
}