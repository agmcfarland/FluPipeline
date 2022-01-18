#alex start:installing packages
#library(remotes)
#install.packages('optparse')
#install_version('latticeExtra','0.6-28') #2016 version
#source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite('ShortRead', suppressUpdates=TRUE, ask=FALSE)
#biocLite('genbankr', suppressUpdates=TRUE, ask=FALSE)
# alex end:installing packages

library(ShortRead)
library(tidyverse)
library(optparse)
library(genbankr)
library(GenomicRanges)

# planet_paragon processing settings.
# Conda enviromnent named fgbio with fgbio installed required.
#planet_paragon_primerTrimmingCoords <- '/home/common/SARS-CoV-2-Philadelphia/data/references/NC_045512.2.primerTrimCoords' #alex comment:code hashed out by alex
#planet_paragon_fgbioCommand   <- 'java -jar /home/everett/miniconda3/envs/fgbio/share/fgbio/fgbio.jar' #alex comment:code hashed out by alex
#planet_paragon_refGenomeFasta <- '/home/common/SARS-CoV-2-Philadelphia/data/references/NC_045512.2.fasta' #alex comment:code hashed out by alex
#planet_paragon_refGenomeBWA   <- '/home/common/SARS-CoV-2-Philadelphia/data/references/NC_045512.2.fasta' #alex comment:code hashed out by alex
#planet_paragon_minAmpliconLength <- 100 #alex comment:code hashed out by alex
#planet_paragon_maxAmpliconLength <- 200 #alex comment:code hashed out by alex

# All options default to the default process_method.
option_list = list(
  make_option(c("--softwareDir"), type="character", default='/home/common/SARS-CoV-2-Philadelphia', help="path to work directory", metavar="character"),
  make_option(c("--process_method"), type="character", default='bushman_artic_v2', help="process method", metavar="character"),
  make_option(c("--workDir"), type="character", default='/tmp', help="path to work directory", metavar="character"),
  make_option(c("--outputFile"), type="character", default='save.RData', help="path to output RData file", metavar="character"),
  make_option(c("--R1"), type="character", default=NULL, help="comma delimited list of R1 fastq files", metavar="character"),
  make_option(c("--R2"), type="character", default=NULL, help="comma delimited list of R2 fastq files", metavar="character"),
  make_option(c("--minAmpliconLength"), type="integer", default=50, help="comma delimited list of R2 fastq files", metavar="character"),
  make_option(c("--maxAmpliconLength"), type="integer", default=350, help="comma delimited list of R2 fastq files", metavar="character"),
  make_option(c("--refGenomeFasta"), type="character", default='/home/common/SARS-CoV-2-Philadelphia/data/references/Wuhan-Hu-1.fasta', help="reference genome FASTA path", metavar="character"),   
  make_option(c("--refGenomeBWA"), type="character", default='/home/common/SARS-CoV-2-Philadelphia/data/references/Wuhan-Hu-1.fasta', help="path to ref genome BWA database", metavar="character"),   
  make_option(c("--refGenomeGenBank"), type="character", default='/home/common/SARS-CoV-2-Philadelphia/data/references/Wuhan-Hu-1.gb', help="path to ref genome genBank file", metavar="character"), 
  make_option(c("--minVariantPhredScore"), type="integer", default=20, help="minimum PHRED score allowed for called varinats", metavar="character"),
  make_option(c("--bwaPath"), type="character", default='/home/everett/ext/bwa', help="path to bwa binary", metavar="character"), 
  make_option(c("--megahitPath"), type="character", default='/home/everett/ext/megahit/bin/megahit', help="path to megahit binary", metavar="character"), 
  make_option(c("--minBWAmappingScore"), type="integer", default=30, help="minimum BWA mapping score", metavar="character"), 
  make_option(c("--minPangolinConf"), type="numeric", default=0.9, help="minimum pangolin confidence value (0-1)", metavar="character"), 
  make_option(c("--samtoolsBin"), type="character", default='/home/everett/ext/samtools/bin', help="path to samtools bin", metavar="character"), 
  make_option(c("--condaShellPath"), type="character", default='/home/everett/miniconda3/etc/profile.d/conda.sh', help="path to conda.sh", metavar="character"),
  make_option(c("--bcftoolsBin"), type="character", default='/home/everett/ext/bcftools/bin',  help="path to bcftools bin", metavar="character"),
  make_option(c("--trimBamCommand"), type="character", default='/home/everett/ext/bamUtil_1.0.15/bam trimBam',  help="Trim bam command", metavar="character"),
  make_option(c("--removeClippedReads"), type="logical", default=FALSE,  help="Number of NTs to remove from the end of alignments", metavar="character"),
  make_option(c("--removeNTsFromAlignmentEnds"), type="integer", default=3,  help="Number of NTs to remove from the end of alignments", metavar="character"),
  make_option(c("--trimQualCode"), type="character", default='5',  help="Min qual trim code", metavar="character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#alex start:FOR TESTING MANUALLY assign opt vals
#opt$softwareDir <- '/home/agmcfarland/flu_project/flu_pipeline'
#opt$process_method <- 'bushman_artic_v2'
#opt$outputFile <- '/home/agmcfarland/flu_project/test/test3/sampleOutputs/IBV_Yamagata_Ref_snpindel/IBV_Yamagata_Ref_snpindel.Rdata'
#opt$workDir <- '/home/agmcfarland/flu_project/test/test3/sampleOutputs/IBV_Yamagata_Ref_snpindel'
#opt$R1 <- '/home/agmcfarland/flu_project/test/test_data/IBV_Yamagata_Ref_snpindel_R1_001.fastq.gz'
#opt$R2 <- '/home/agmcfarland/flu_project/test/test_data/IBV_Yamagata_Ref_snpindel_R2_001.fastq.gz'
#opt$refGenomeBWA <- '/home/agmcfarland/flu_project/test/test_data/IBV_Yamagata_Ref.fasta'
#opt$refGenomeFasta <- '/home/agmcfarland/flu_project/test/test_data/IBV_Yamagata_Ref.fasta' #alex comment:same as bwa.
#opt$refGenomeGenBank <- '/home/agmcfarland/flu_project/test/test_data/IBV_Yamagata_Ref.gb'
#alex end:FOR TESTING MANUALLY assign opt vals

#alex start:set working directory to workDir
setwd(opt$workDir)
#alex end:set working directory to workDir

#alex start:copy bwa files to workDir and change the reference path locations to the files in workDir
system(paste0('cp ',opt$refGenomeBWA, ' .'))
system(paste0('cp ',opt$refGenomeFasta, ' .'))
#system(paste0('cp ',opt$refGenomeGenBank, ' .'))
opt$refGenomeBWA <- paste0(opt$workDir,'/',basename(opt$refGenomeBWA))
opt$refGenomeFasta <- paste0(opt$workDir,'/',basename(opt$refGenomeFasta))
#opt$refGenomeGenBank <- paste0(opt$workDir,'/',basename(opt$refGenomeGenBank))
#alex end:copy bwa files to workDir and change the reference path locations to the files in workDir

#alex start:MANUAL change source
source(paste0(opt$softwareDir, '/lib/assemble.lib.R'))#source(paste0(opt$softwareDir, '/lib/assemble.lib.R'))
#alex end:change source


# Set up non-default processing parameter
if(opt$process_method == 'planet_paragon'){
  opt$removeNTsFromAlignmentEnds <- 0
  opt$removeClippedReads <- FALSE
  opt$refGenomeFasta <- planet_paragon_refGenomeFasta
  opt$refGenomeBWA   <- planet_paragon_refGenomeBWA
  opt$fgbioCommand   <- planet_paragon_fgbioCommand
  opt$minAmpliconLength <- planet_paragon_minAmpliconLength
  opt$maxAmpliconLength <- planet_paragon_maxAmpliconLength
  opt$fgbioTrimPrimerCoords <- planet_paragon_primerTrimmingCoords
}

if(opt$process_method == 'bushman_artic_v2'){
  opt$removeNTsFromAlignmentEnds <- 3
  opt$removeClippedReads <- FALSE #alex comment:changed this from TRUE to FALSE
}


# Handle missing required parameters.
if(! 'R1' %in% names(opt)) stop('--R1 must be defined.')
if(! 'R2' %in% names(opt)) stop('--R2 must be defined.')


# Create work directory.
#if(dir.exists(opt$workDir)) stop('Error -- work directory already exists') #alex comment:hashed out by alex
#dir.create(opt$workDir) #alex comment:hashed out by alex
#if(! dir.exists(opt$workDir)) stop('Error -- could not create the work directory.') #alex comment:hashed out by alex

# Create table of software and R package version numbers.
opt$softwareVersionTable <- createSoftwareVersionTable()


# Define data structures which must be present even if empty for downstream analyses.
opt$errorCode <- 0
opt$errorMessage <- NA
opt$pangolinAssignment <- NA
opt$pangolinAssignmentConflict <- NA
opt$pangolinAssignmentPangoLEARN_version <- NA
opt$variantTable      <- data.frame()
opt$variantTableMajor <- data.frame()
opt$contigs           <- Biostrings::DNAStringSet()


R1s <- unlist(strsplit(opt$R1, ','));  if(! all(file.exists(R1s))) stop('All the R1 files could not be found.')
R2s <- unlist(strsplit(opt$R2, ','));  if(! all(file.exists(R2s))) stop('All the R1 files could not be found.')

t1 <- paste0(opt$workDir, '/x')

# Combine R1 and R2 data files in composite R1 and R2 files.
system(paste0('cat ', paste0(R1s, collapse = ' '), ' > ', t1, '_R1.fastq'))
system(paste0('cat ', paste0(R2s, collapse = ' '), ' > ', t1, '_R2.fastq'))


# Quality trim reads and create trimmed FASTA files.
r <- prepareTrimmedReads(readFastq(paste0(t1, '_R1.fastq')), readFastq(paste0(t1, '_R2.fastq')), qualCode = opt$trimQualCode)
writeFasta(r[[1]], file = paste0(t1, '_R1.trimmed.fasta'))
writeFasta(r[[2]], file = paste0(t1, '_R2.trimmed.fasta'))


# Align trimmed reads to the reference genome.
#alex start:make bwa index file for reference
system(paste0(opt$bwaPath, ' index ', opt$refGenomeBWA)) #makes index files using the whole filename as the prefix
#alex end:make bwa index file for reference
system(paste0(opt$bwaPath, ' mem -M ', opt$refGenomeBWA, ' ',  paste0(t1, '_R1.trimmed.fasta'), ' ', 
              paste0(t1, '_R2.trimmed.fasta'),  ' > ', paste0(t1, '_genome.sam')))

# Prevent the alignment of clipped reads.
# Here we align the quality trimmed reads to the reference genome, identify which reads are clipped, 
# create new FASTA inputs (without reads that were clipped), and then align again. Running bwa twice
# is not ideal, a better solution should be developed.
#----------------------------------------------------------------------------------------------------
if(opt$removeClippedReads){
  message('Removing clipped reads\n') #alex comment:this step takes a long time, mostly due to variable `o` being initialized.
  
  system(paste0("grep -v '^@' ", paste0(t1, '_genome.sam'), "> ",paste0(t1, '_genome.samNoHeader')))
  o <- read.table(paste0(t1, '_genome.samNoHeader'), sep = '\t', header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
  
  a <- o[! grepl('H|S', toupper(o$V6)),] #alex comment:more info on H,S designations in CIGAR string https://github.com/GregoryFaust/samblaster/issues/35
  b <- o[grepl('H|S', toupper(o$V6)),]
  write(unique(b$V1), file = paste0(t1, '_alignmentClippedReads')) #alex comment:all reads that were clipped
  
  opt$percentClippedReadsPairsRemoved <- sprintf("%.2f%%", (n_distinct(b$V1)/n_distinct(o$V1))*100)
  
  R1 <- readFasta(paste0(t1, '_R1.trimmed.fasta'))
  R2 <- readFasta(paste0(t1, '_R2.trimmed.fasta'))

  writeFasta(R1[! as.character(R1@id) %in% unique(b$V1)], file = paste0(t1, '_R1.trimmed.fasta'), mode = 'w')#alex comment:reads that are not clipped are overwritten
  writeFasta(R2[! as.character(R2@id) %in% unique(b$V1)], file = paste0(t1, '_R2.trimmed.fasta'), mode = 'w')
  
  system(paste0(opt$bwaPath, ' mem -M ', opt$refGenomeBWA, ' ',  paste0(t1, '_R1.trimmed.fasta'), ' ', 
                paste0(t1, '_R2.trimmed.fasta'),  ' > ', paste0(t1, '_genome.sam'))) #alex comment:second alignment step
}


system(paste0(opt$samtoolsBin, '/samtools view -S -b ', paste0(t1, '_genome.sam'), ' > ', paste0(t1, '_genome.bam'))) #alex comment:sam to bam conversion
invisible(file.remove(paste0(t1, '_genome.sam')))

# Remove non-paired mates and filter alignment by insert size.

comm <- paste0(opt$samtoolsBin, '/samtools view -F 8 -h ', t1, "_genome.bam | awk 'substr($0,1,1)==\"@\" || ($9 >= ", opt$minAmpliconLength,
               " && $9 <= ", opt$maxAmpliconLength,  ") || ($9 <= -", opt$minAmpliconLength, " && $9 >=- ", opt$maxAmpliconLength, ")' | ",
               opt$samtoolsBin, '/samtools view -b > ', t1, '_genome.filt.bam') #alex comment:come back to this and look more closely. a lot information here.
system(comm)


#system(paste0(' cp ',  t1, '_genome.bam ',  t1, '_genome.filt.bam')) #alex comment:this line was hashed out in the original script

#alex start:commented out planet paragon block
#if(opt$process_method == 'planet_paragon'){
#comm <- paste0(opt$fgbioCommand, ' TrimPrimers -i ', t1, '_genome.filt.bam -o ', t1, '_genome.filt.fgbio.bam -H true -p ', opt$fgbioTrimPrimerCoords)

#p <- c('#!/bin/bash', paste0('source ', opt$condaShellPath), 'conda activate fgbio', comm)
#writeLines(p, file.path(opt$workDir, 'fgbio.script'))
#system(paste0('chmod 755 ', file.path(opt$workDir, 'fgbio.script')))
#system(file.path(opt$workDir, 'fgbio.script'))

# Copy the primer trimmed to the orignal file so the pipeline can continue.
#system(paste0('cp ', t1, '_genome.filt.bam ',  t1, '_genome.filt.org.bam'))
#system(paste0('cp ', t1, '_genome.filt.fgbio.bam ', t1, '_genome.filt.bam'))
#opt$fgbioPrimerTrimmed <- TRUE
#}
#alex end:commented out planet paragon block

# Trim NTs from the ends of all aligned reads.#alex comment:this step is suggested by John to improve SNP detection, espcially with deletions
#alex comment: it doesn't seem like there is an option to not do this either.
#-------------------------------------------------------------------------------------------------- 
if(opt$removeNTsFromAlignmentEnds > 0){
  message('Trimming alignments by ', opt$removeNTsFromAlignmentEnds > 0, ' NTs\n')
  system(paste0(opt$trimBamCommand, ' ', t1, '_genome.filt.bam ', t1, '_genome.filt.endTrim.bam ', opt$removeNTsFromAlignmentEnds, ' --clip'))
  system(paste0('cp ', t1, '_genome.filt.endTrim.bam ', t1, '_genome.filt.bam'))
}


# Remove read pairs with mapping qualities below the provided min score.
system(paste0(opt$samtoolsBin, '/samtools view -q ', opt$minBWAmappingScore, ' -b ', t1, '_genome.filt.bam > ', t1, '_genome.filt.qual.bam'))


# Retrieve a list of aligned reads.
alignedReadsIDs <- system(paste0(opt$samtoolsBin, '/samtools view ', t1, '_genome.filt.qual.bam | cut  -f 1 | uniq'), intern = TRUE)

if(length(alignedReadsIDs) == 0){
  opt$errorCode <- 1
  opt$errorMessage <- 'No quality trimmed reads aligned to the reference genome.'
  save(opt, file = opt$outputFile)
  #unlink(opt$workDir, recursive = TRUE) #alex comment:hashed out by alex.
  stop()
}

#alex start:skipping denovo assembly step
#code would be here
#alex end: skipping denovo assembly step


# Sort and index filtered genome alignment.
system(paste0(opt$samtoolsBin, '/samtools sort -o ', paste0(t1, '_genome.filt.qual.sorted.bam'), ' ', paste0(t1, '_genome.filt.qual.bam')))
system(paste0(opt$samtoolsBin, '/samtools index ', paste0(t1, '_genome.filt.qual.sorted.bam')))


# Determine the maximum read depth. Overlapping mates will count a position twice.
system(paste0(opt$samtoolsBin, '/samtools depth -d 0 ', paste0(t1, '_genome.filt.qual.sorted.bam'), ' -o ', paste0(t1, '.depth'))) #alex comment:depth file shows all segments have per base depth reported
maxReadDepth <- max(read.table(paste0(t1, '.depth'), sep = '\t', header = FALSE)[,3])


# Create pileup data file for determining depth at specific locations.
# (!) mpileup will remove duplicate reads.
#alex comment:http://www.htslib.org/doc/samtools-mpileup.html
#alex comment:-A(count orphans) -a(output all positions, even those with 0 depth) -Q(minimum base quality for a base to be considered) -d(max read depth to filter at)
system(paste0(opt$samtoolsBin, '/samtools mpileup -A -a -Q 0 -o ', paste0(t1, '.pileup'), ' -d ', maxReadDepth, 
              ' -f ', opt$refGenomeFasta, ' ', paste0(t1, '_genome.filt.qual.sorted.bam')))


# Determine the percentage of the reference genome covered in the pileup data.
#alex comment:don't try to print this variable in rstudio, it is huge
opt$pileupData <- tryCatch({
  read.table(paste0(t1, '.pileup'), sep = '\t', header = FALSE, quote = '')[,1:5]
}, error = function(e) {
  return(data.frame())
})



# Pileup format reports the number of read pairs (column 4) while VCF format (DP) 
# reports the number of reads which appears to report 2x the pileup format value. 
# Confirmed by looking at pileup in IGV.


if(nrow(opt$pileupData) > 0){
  #alex start:reconfiguring pileup coverge calculation to work with multiple segments
  #refGenomeLength <- nchar(as.character(readFasta(opt$refGenomeFasta)@sread)) #original code hashed out by alex
  #opt$refGenomePercentCovered <- nrow(subset(opt$pileupData,  V4 >= 1))  / refGenomeLength #original code hashed out by alex
  #opt$refGenomePercentCovered_5reads <- nrow(subset(opt$pileupData,  V4 >= 5))  / refGenomeLength #original code hashed out by alex
  #start new
  df_sample <- opt$pileupData%>%
    dplyr::select(V1,V2,V4)%>%#select segment,baseposition,read coverage
    dplyr::group_by(V1)%>%
    dplyr::summarize(ref_genome_length = n(),
                     gt_1=sum(V4 >= 1),
                     gt_5=sum(V4 >= 5)
    )%>%
    dplyr::mutate(perc_gt_1 = gt_1/ref_genome_length*100,
                  perc_gt_5 = gt_5/ref_genome_length*100)
  opt$refGenomePercentCovered <- df_sample%>%select(V1,perc_gt_1)
  opt$refGenomePercentCovered_5reads <- df_sample%>%select(V1,perc_gt_5)
  # these variables are now tibbles and not single numeric values
  #alex end:reconfiguring pileup coverge calculation to work with multiple segments
  
  #alex comment:variant calling in the line below. mpileup command piped into call command.
  system(paste0(opt$bcftoolsBin, '/bcftools mpileup -A -Ou -f ', opt$refGenomeFasta, ' -d 10000 -L 10000 ',
                paste0(t1, '_genome.filt.qual.sorted.bam'), ' |  ', opt$bcftoolsBin,  '/bcftools call -mv -Oz ', 
                ' -o ', paste0(t1, '.vcf.gz')))
  
  #opt$minVariantPhredScore <- 2 #alex comment:added by alex for testing purposes
  # Read in the variant table created by bcf tools. 
  # We use tryCatch() here because the table may be empty only containing the header information.
  #alex comment:DP(Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases) that support that variant call
  opt$variantTable <- tryCatch({
    system(paste0(opt$bcftoolsBin, "/bcftools filter -i'QUAL>", opt$minVariantPhredScore, " && DP>4' ", 
                  paste0(t1, '.vcf.gz'), " -O z -o ", paste0(t1, '.filt.vcf.gz')))
    
    system(paste0(opt$bcftoolsBin, '/bcftools index ', t1, '.filt.vcf.gz'))
    
    x <- read.table(paste0(t1, '.filt.vcf.gz'), sep = '\t', header = FALSE, comment.char = '#')
    names(x) <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'OTHER')
    x <- x[, !is.na(names(x))]
    x #alex comment:returns x as opt$variantTable
  },  error=function(cond) {
    return(data.frame()) 
  })
  
  
  #x <- read.table(paste0(t1, '.vcf.gz'), sep = '\t', header = FALSE, comment.char = '#') #added by alex for testing
  #names(x) <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'OTHER')#added by alex for testing
  #x <- x[, !is.na(names(x))]#added by alex for testing
  
  if(nrow(opt$variantTable) > 0){
    
    opt$variantTable <- tryCatch({
      # Here we parse the pileup data to create a more informative alt call for variants.
      #alex comment:for each row in the variant table generated by `bcftools call` and reference that row in the pileup data generated by `samtools pileup`.
      #alex comment:In the pileup data extract the string from the column containing the per read base stats http://www.htslib.org/doc/samtools-mpileup.html (keyword: This encodes information matches,)
      #alex comment:parsePileUpString2 counts the number of base substitutions observed uses grepl or indel info using a boolean flag to count indel length and sequence. It's a little complicated.
      x <- bind_rows(lapply(1:nrow(opt$variantTable), function(i){
        x <- opt$variantTable[i,]
        
        ###if(x$POS == 25445) browser()
        #d <- subset(opt$pileupData, V2 == x$POS) #original code hashed out by alex
        #alex start:adding line for multi chromosome support
        d <- subset(opt$pileupData, V2 == x$POS & V1 == as.character(x$CHROM))
        #alex end:adding line for multi chromosome support
        p <- parsePileUpString2(as.character(d$V5))
        
        # Expand the variant call to include the different possibilities.
        x <- x[rep(1, length(p)),] #alex comment: make new identical rows of x (which should be one row) that are the length of p. Then below alter it to contain additional variant information.
        x$ALT <- names(p)
        x$altCount <- as.numeric(p)
        
        dplyr::select(x, -ID, -FILTER, -INFO, -FORMAT, -OTHER) #alex comment: left with these columns: CHROM, POS, REF, ALT, QUAL, altCount
      }))
      x #alex comment:returns x as opt$variantTable
    }, error=function(cond) {
      #alex start:adding error message to this stop trycatch
      opt$errorCode <- 5
      opt$errorMessage <- 'Error parsing variant occurrences.'
      save(opt, file = opt$outputFile)
      #alex end:adding error message to this stop trycatch
      stop(paste0(opt$workDir,' Error parsing variant occurrences.'))
    })
    
    # Increment the position of indels because their calls are from the NT preceding the event.
    opt$variantTable$POS <- ifelse(grepl('ins|del', opt$variantTable$ALT), opt$variantTable$POS+1, opt$variantTable$POS)
    
    # Retrieve position read counts from pileup data. This is done here rather than previously because we 
    # had to increment the indel positions.
    
    #alex start:make the process of merging updated read base position in opt$variantTable with read count pileup data in opt$pileupData able to handle multiple chromosomes
    
    #opt$variantTable <- bind_rows(lapply(split(opt$variantTable, opt$variantTable$POS), function(x){ #original code hashed out by alex
    #browser() #original code hashed out by alex
    #x$reads <- subset(opt$pileupData, V2 == x$POS[1])$V4 #original code hashed out by alex
    #x$percentAlt <- x$altCount / x$reads[1] #original code hashed out by alex
    #x #original code hashed out by alex
    #})) #original code hashed out by alex
    
    #start new
    # important: the rownames are different, but i don't think this should affect anything downstream
    opt$variantTable <- merge(opt$variantTable, opt$pileupData%>%select(-V3,-V5), 
                              by.x=c('CHROM','POS'), 
                              by.y=c('V1','V2'), #alex comment:this is the same as CHROM and POS in opt$variantTable
                              all.x=TRUE) #alex comment:enforce that nrow is the same before and after merger
    opt$variantTable <- opt$variantTable%>%
      dplyr::rename(reads=V4)%>%
      dplyr::mutate(percentAlt=altCount/reads)
    #alex end:merge updated read base position in opt$variantTable with read count pileup data in opt$pileupData 
    
    
    # Here we undo the correction of the indel positions because the following code expects uncorrected positions.
    # The positions will be re-corrected at the end.
    opt$variantTable$POS <- ifelse(grepl('ins|del', opt$variantTable$ALT), opt$variantTable$POS-1, opt$variantTable$POS)
    
    
    # Select the major variant for each position.
    opt$variantTableMajor <- dplyr::filter(opt$variantTable,  percentAlt > 0.5 & ! ALT == 'e') %>%
      #dplyr::group_by(POS) #original code hashed out by alex
      dplyr::group_by(CHROM,POS) %>% #alex comment: updated this line so that calculation takes multiple chromosomes into account
      dplyr::top_n(1, wt = percentAlt) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::filter(reads > 4)
    
  }
} else {
  opt$refGenomePercentCovered <- 0
  opt$refGenomePercentCovered_5reads  <- 0
  
  opt$errorCode <- 5
  opt$errorMessage <- 'No pileup or variant data available.'
}

#alex start:hashed out this entire code block for predicting effect of SNPs/INDELs on viral proteins

##### SHIFTING THE DEL VALUES ARE PREVENTING A MATCH IN THE VCF FILES.


# # Determine the result of the variants on the AA sequence of viral proteins.
# # This approach assumes all orfs are non-overlapping pos. strand sequenes.
# if(nrow(opt$variantTableMajor) > 0){

#   # Here we copy the variant vcf file and selectively remove calls 
#   # which are not found in our filtered variantTable. This is done 
#   # to preserve the additional comment lines which appear to be necessary. 

#   tryCatch({
#     system(paste0('cp ', t1, '.filt.vcf.gz ', t1, '.filt.vcf.copy.gz'))
#     system(paste0('gunzip ', t1, '.filt.vcf.copy.gz'))

#     o <- lapply(readLines(paste0(t1, '.filt.vcf.copy')), function(x){
#       if(grepl('^#', x)){
#         return(x)
#       } else {
#         if(as.integer(unlist(strsplit(x, '\t'))[2]) %in% opt$variantTableMajor$POS){
#           return(x)
#         } else {
#           return(NULL)
#         }
#       }
#     })

#     write(unlist(o[! sapply(o, is.null)]), file = paste0(t1, '.filt.vcf.copy2'))

#     # Capture vcf
#     # read.table(textConnection(unlist(opt$finalVCF)), sep = '\t')
#     opt$finalVCF = o[! sapply(o, is.null)]

#     system(paste0('bgzip ', t1, '.filt.vcf.copy2'))
#     system(paste0(opt$bcftoolsBin, '/bcftools index ', t1, '.filt.vcf.copy2.gz'))
#     system(paste0('cat  ', opt$refGenomeFasta, ' | ', opt$bcftoolsBin, '/bcftools consensus ',  
#                   t1, '.filt.vcf.copy2.gz > ', t1, '.consensus.fasta'))

#     opt$concensusSeq <- as.character(readFasta(paste0(t1, '.consensus.fasta'))@sread)

#     # Predict pangolin lineage for genomes with >= 90% coverage.
#     if(opt$refGenomePercentCovered_5reads >= 0.90){

#       # Create a bash script which will start the requires Conda environment and run pangolin.
#       p <- c('#!/bin/bash', paste0('source ', opt$condaShellPath), 'conda activate pangolin', 'pangolin -o $2 $1')
#       writeLines(p, file.path(opt$workDir, 'pangolin.script'))
#       system(paste0('chmod 755 ', file.path(opt$workDir, 'pangolin.script')))

#       comm <- paste0(file.path(opt$workDir, 'pangolin.script'), ' ', t1, '.consensus.fasta ', t1, '.consensus.pangolin')
#       system(comm)

#       if(file.exists(paste0(t1, '.consensus.pangolin/lineage_report.csv'))){
#         o <- read.csv(paste0(t1, '.consensus.pangolin/lineage_report.csv')) 

#         if(nrow(o) > 0){
#           opt$pangolinAssignment <- as.character(o[1,]$lineage)
#           opt$pangolinAssignmentConflict <- as.numeric(o[1,]$conflict)
#           opt$pangolinAssignmentPangoLEARN_version <- as.character(o[1,]$pangoLEARN_version)
#         }
#       }
#     } 
#   }, error=function(cond) {
#     stop('Error creating concensus sequence.')
#   })

#   message('Concensus sequence is ', refGenomeLength - nchar(opt$concensusSeq), ' NT shorter than the reference sequence.')

#   gb <- readGenBank(opt$refGenomeGenBank)

#   cds <- gb@cds
#   seqlevels(cds) <- 'genome'
#   seqnames(cds)  <- 'genome'

#   # Calculate how the shift left or right caused by deletions and insertions.
#   opt$variantTableMajor$shift <- ifelse(grepl('del', opt$variantTableMajor$ALT), (nchar(opt$variantTableMajor$ALT)-3)*-1, 0)
#   opt$variantTableMajor$shift <- ifelse(grepl('ins', opt$variantTableMajor$ALT), (nchar(opt$variantTableMajor$ALT)-3), opt$variantTableMajor$shift)


#   # Remove variant positions flanking indels since they appear to be artifacts. 
#   artifacts <- c(opt$variantTableMajor[grep('ins|del', opt$variantTableMajor$ALT),]$POS + abs(opt$variantTableMajor[grep('ins|del', opt$variantTableMajor$ALT),]$shift),
#                  opt$variantTableMajor[grep('ins|del', opt$variantTableMajor$ALT),]$POS -1)

#   opt$variantTableMajor <- opt$variantTableMajor[! opt$variantTableMajor$POS %in% artifacts,]

#   opt$variantTableMajor <- bind_rows(lapply(split(opt$variantTableMajor, 1:nrow(opt$variantTableMajor)), function(x){

#     # Determine the offset of this position in the concensus sequence because it may not be the same length
#     # if indels have been applied. Here we sum the indel shifts before this variant call.

#     offset <- sum(opt$variantTableMajor[1:grep(x$POS, opt$variantTableMajor$POS),]$shift) 

#     cds2 <- cds
#     start(cds2) <- start(cds2) + offset 
#     end(cds2) <- end(cds2) + offset 

#     v1 <- GRanges(seqnames = 'genome', ranges = IRanges(x$POS, end = x$POS), strand = '+')
#     o1 <- GenomicRanges::findOverlaps(v1, cds)

#     v2 <- GRanges(seqnames = 'genome', ranges = IRanges(x$POS + offset, end = x$POS + offset), strand = '+')
#     o2 <- GenomicRanges::findOverlaps(v2, cds2)

#     if(length(o2) == 0){
#       x$genes <- 'intergenic'

#       if (grepl('ins', as.character(x$ALT))){
#         x$type <- paste0('ins ', nchar(x$ALT)-3)
#       } else if (grepl('del', as.character(x$ALT))){
#         x$type <- paste0('del ', nchar(x$ALT)-3)
#       } else {
#         x$type <- ' '
#       }
#     } else {

#       # Define the gene the variant is within.
#       hit1 <- cds[subjectHits(o1)]
#       hit2 <- cds2[subjectHits(o2)]

#       x$genes <- paste0(hit2$gene, collapse = ', ')

#       # Native gene AA sequence.
#       orf1  <- as.character(translate(DNAString(substr(as.character(readFasta(opt$refGenomeFasta)@sread), start(hit1), end(hit1)))))

#       # Variant gene AA sequence.
#       orf2 <- as.character(translate(DNAString(substr(opt$concensusSeq, start(hit2), end(hit2)))))


#       # Determine the offset of this position in the concensus sequence because it may not be the same length
#       # if indels have been applied. Here we sum the indel shifts before this variant call.
#       # offset <- sum(opt$variantTableMajor[1:grep(x$POS, opt$variantTableMajor$POS),]$shift)

#       #              1   2   3   4   5   6   7   8
#       # 123 456 789 012 345 678 901 234 567 890 123
#       # ATG CAT TGA ATG GGC TTA CGA GCT TAA GTA TAG
#       #             ^             x  21-10 + 2 = 13/3 = 4.3 ~ 4
#       #                          x   20-10 + 2 = 12/3 = 4.0 = 4
#       #                         x    19-10 + 2 = 11/3 = 3.6 ~ 4
#       #                                 x   25-10 + 2 = 17/3 = 5.6 ~ 6
#       #                                  x  26-10 + 2 = 18/3 = 6.0 = 6
#       #                                   x 27-10 + 2 = 19/3 = 6.3 ~ 6

#       aa <- round(((x$POS - start(hit1)) + 2)/3)
#       orf_aa <- substr(orf1, aa, aa)

#       aa2 <- round((((x$POS + offset) - start(hit2)) + 2)/3)
#       orf2_aa <- substr(orf2, aa2, aa2)

#       maxALTchars <- max(nchar(unlist(strsplit(as.character(x$ALT), ','))))

#       if(nchar(as.character(x$REF)) == 1 & nchar(as.character(x$ALT)) > 1 & maxALTchars == 1){
#         x$type <- paste0(x$POS, '_mixedPop')
#       } else if (grepl('ins', as.character(x$ALT))){
#         x$type <- paste0('ins ', nchar(x$ALT)-3)
#       } else if (grepl('del', as.character(x$ALT))){
#         x$type <- paste0('del ', nchar(x$ALT)-3)
#       } else if (orf_aa != orf2_aa){
#         # JKE 2021-08-18.
#         # Fixes probject of multiple hits.
#         ###x$type <- paste0(orf_aa, aa2, orf2_aa)
#         x$type <- paste0(orf_aa[1], aa2[1], orf2_aa[1])
#       } else {
#         x$type <- 'silent'
#       }
#     }
#     x
#   }))

# } else {
#   # There were no variants called so we report the reference as the concensus sequence.
#   opt$concensusSeq <- as.character(readFasta(opt$refGenomeFasta)@sread)
# }

#alex end:hashed out this entire code block for predicting effect of SNPs/INDELs on viral proteins


# BCFtools calls indels by the base preceding the modification.
# Here we find deletion calls and increment the position by one to mark the first deleted base.
i <- opt$variantTable$POS %in% opt$variantTable[grep('del', opt$variantTable$ALT),]$POS
if(any(i)) opt$variantTable[i,]$POS <- opt$variantTable[i,]$POS + 1

i <- opt$variantTableMajor$POS %in% opt$variantTableMajor[grep('del', opt$variantTableMajor$ALT),]$POS
if(any(i)) opt$variantTableMajor[i,]$POS <- opt$variantTableMajor[i,]$POS + 1

write.csv(opt$variantTable,'variantTable.csv',row.names=FALSE)
write.csv(opt$variantTableMajor,'variantTableMajor.csv',row.names=FALSE)

save(opt, file = opt$outputFile)
#unlink(opt$workDir, recursive = TRUE) #alex comment:delete all files in the folder


