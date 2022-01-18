
createGisadUploadForm <- function(genomeMetaData, gisaideUser){
  if(file.exists('summaries/genomes/GISAD_upload.xlsx')) file.remove('summaries/genomes/GISAD_upload.xlsx')
  openxlsx::write.xlsx(
    data.frame(Submitter = rep(gisaideUser, nrow(genomeMetaData)),
               'FASTA filename' = 'x.fasta',
               'Virus name' = genomeMetaData$genome_id,
               'Type' = 'betacoronavirus',
               'Passage details/history' = 'original',
               'Collection date' = genomeMetaData$sample_date,
               'Location' = paste0('North America/USA/', genomeMetaData$state),
               'Additional location information' = NA,
               'Host' = 'Human',
               'Additional host information' = NA,
               'Sampling Strategy' = genomeMetaData$covv_sampling_strategy,
               'Gender' = 'unknown',
               'Patient age' = 'unknown',
               'Patient status' = 'unknown',
               'Specimen source' = NA,
               'Outbreak' = NA,
               'Last vaccinated' = NA,
               'Treatment' = NA,
               'Sequencing technology' = 'Illumina NextSeq', 
               'Assembly method' = 'Samtools v. 1.10; BCFtools v. 1.10.2',
               'Coverage' = paste0(genomeMetaData$mean_coverage, 'x'),
               'Originating lab' = 'Hospital of the University of Pennsylvania Molecular Pathology Lab',
               'Address' = '3400 Spruce St. Philadelphia, PA 19104',
               'Sample ID given by originating laboratory' = genomeMetaData$lab_id,
               'Submitting lab' = 'Bushman Lab - University of Pennsylvania', 
               'Address' = '425 Johnson Pavilion 3610 Hamilton Walk Philadelphia, PA 19104',
               'Sample ID given by the submitting laboratory' = NA,
               'Authors' =  'John K. Everett, Kyle Rodino, Andrew D. Marques, Shantan Reddy, Aoife M. Roche, Young Hwang, Scott Sherrill-Mix, Samantha A. Whiteside, Jevon Graham-Wooten, Layla A. Khatib, Ayannah S. Fitzgerald, Arupa Ganguly, Mike Feldman, Brendan Kelly, Ronald G. Collman and Frederic Bushan',
               'Comment' = NA,
               'Comment Icon' = NA, check.names = FALSE), file = 'summaries/genomes/GISAD_upload.xlsx')
}


createNIHUploadForm <- function(genomeMetaData){
  if(file.exists('summaries/genomes/NIH_upload.xlsx')) file.remove('summaries/genomes/NIH_upload.xlsx')
  openxlsx::write.xlsx(
    data.frame(Sequence_ID = genomeMetaData$genome_id,
               country = paste0('USA: ', genomeMetaData$stateName, ', ', genomeMetaData$city),
               host = 'Homo sapiens',
               isolate = genomeMetaData$lab_id,
               'collection-date' = genomeMetaData$sample_date, check.names = FALSE), file = 'summaries/genomes/NIH_upload.xlsx')
}


tmpFile <- function(){ paste0('tmp.', paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = '')) }

removeProblematicSamples <- function(x){
  if('exp' %in% names(x))        x <- dplyr::filter(x, ! grepl('VSP0069', exp))
  if('patient_id' %in% names(x)) x <- dplyr::filter(x, patient_id != '237')
  
  if('trial_id' %in% names(x)) x <- subset(x, ! trial_id %in% c('CHOP_Planet',  'Lennnon_animals', 'simulated_dataset'))
  
  # Use Testing date rather than collection date where available.
  if('sampleCollection_date' %in% names(x) & 'Testing_Date' %in% names(x)){
    x$sampleCollection_date <- mdy(x$sampleCollection_date)
    x$Testing_Date <- mdy(x$Testing_Date)
    x[! is.na(x$Testing_Date),]$sampleCollection_date <- x[! is.na(x$Testing_Date),]$Testing_Date
  }
  
  x
}


variantSchematic <- function(VSPobjPath){
  library(ggplot2)
  library(tibble)
  
  load(VSPobjPath)
  d <- data.frame(x1 = 0, x2 = nchar(opt$concensusSeq), y1 = 0, y2 = 1)
  m <- data.frame(x1 = opt$variantTableMajor$POS, x2 = opt$variantTableMajor$POS, 
                  y1 = 0, y2 = 1, ALT = as.character(opt$variantTableMajor$ALT))

  m$ALT2 <- ifelse(m$ALT == 'A', 'A', 
                   ifelse(m$ALT == 'T', 'T',
                          ifelse(m$ALT == 'C', 'C',
                                 ifelse(m$ALT == 'G', 'G',
                                        ifelse(grepl('del', m$ALT, ignore.case = T), 'Del',
                                               ifelse(grepl('ins', m$ALT, ignore.case = T), 'Ins', 'Other'))))))

  m$ALT2 <- factor(as.character(m$ALT2), levels = c('A', 'T', 'C', 'G', 'Del', 'Ins', 'Other'))

  p <- ggplot() +
    scale_color_manual(values = c('red', 'blue', 'green2', 'gold2', 'black', 'gray50', 'purple')) +
    geom_rect(data = d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="black", fill = 'white') +
    geom_segment(data = m, aes(x = x1, y = y1, xend = x2, yend = y2, color = ALT2), size = 1) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())

  VSPid <- stringr::str_extract(VSPobjPath, 'VSP\\d+')
  
  ggsave(p, file = file.path('summaries', 'images', paste0(VSPid, '.png')), units = 'in', height = 0.5, width = 3)
  
  o <- subset(samples, VSP == VSPid)
  if(nrow(o) == 1){
    return(tibble('Sample' = VSPid, Date = o$sampleCollection_date, 
                  'Num. mutations' = nrow(opt$variantTableMajor), lineage = as.character(opt$pangolinAssignment), class = o$Rationale,
                  'Analysis report' = paste0('<a href="', file.path('summaries', 'trials', o$trial_id, paste0(o$patient_id, '.pdf')), '">link</a>'),
                  'Mutation positions' = paste0("<img src='",  paste0('summaries/images/', VSPid, '.png'), "'>")))
  } else {
    return(tibble())
  }
}



# Select a sample to represent each subject, sample type, time point combination.
# Select composite samples when available otherwise select the samples with the greated coverage.
# Break ties with read coverage percentages.

representativeSampleSummary <- function(summary, minPercentRefReadCoverage5){
  bind_rows(lapply(split(summary, summary$VSP), function(x){
    top_n(x, 1, wt = percentRefReadCoverage5) %>% dplyr::slice(1)
  })) %>% dplyr::filter(percentRefReadCoverage5 >= minPercentRefReadCoverage5)
}


# Retrieve consensus sequences from representative samples.
retrieveConcensusSeqs <- function(summary){
  cluster <- makeCluster(CPUs)
  clusterExport(cluster, 'softwareDir')

  r <- parLapply(cluster, 1:nrow(summary), function(x){
    library(Biostrings)
    x <- summary[x,]
    f <- paste0(softwareDir, '/summaries/VSPdata/', x$exp, '.', ifelse(x$type == 'composite', 'composite', 'experiment'), '.RData')
    if(! file.exists(f)){
      message('Can not locate VSP data file -- ', f)
      return(DNAStringSet())
    }
    load(f)
    d <- DNAStringSet(opt$concensusSeq)
    
    names(d) <- gsub('\\s+', '_', paste0(x$trial, '|', x$Subject, '|', x$type, '|', 
                                         gsub('\\-', '', x$sampleDate), '|', x$VSP, '|', 
                                         ceiling(mean(opt$pileupData$V4)), 'x mean coverage|', 
                                         as.character(opt$pangolinAssignment)))
    d
  })
  
  stopCluster(cluster)
  Reduce('append', r)
}
