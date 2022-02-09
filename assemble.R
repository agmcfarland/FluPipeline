library(ShortRead)
library(tidyverse)
library(optparse)
library(genbankr)
library(GenomicRanges)

message('start assemble.R')

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
#opt$softwareDir <- '/home/agmcfarland/flu_project/FluPipeline'
#opt$process_method <- 'bushman_artic_v2'
#opt$outputFile <- '/home/agmcfarland/flu_project/FluPipeline/run_test/output/sampleOutputs/H1N1_A_Brisbane_59_2007_snp/H1N1_A_Brisbane_59_2007_snp.Rdata'
#opt$workDir <- '/home/agmcfarland/flu_project/FluPipeline/run_test/output/sampleOutputs/H1N1_A_Brisbane_59_2007_snp'
#opt$R1 <- '/home/agmcfarland/flu_project/FluPipeline/run_test/data/H1N1_A_Brisbane_59_2007_snp_R1_001.fastq.gz'
#opt$R2 <- '/home/agmcfarland/flu_project/FluPipeline/run_test/data/H1N1_A_Brisbane_59_2007_snp_R2_001.fastq.gz'
#opt$refGenomeBWA <- '/home/agmcfarland/flu_project/FluPipeline/run_test/output/sampleOutputs/H1N1_A_Brisbane_59_2007_snp/H1N1_A_Brisbane_59_2007.fasta'
#opt$refGenomeFasta <- '/home/agmcfarland/flu_project/FluPipeline/run_test/output/sampleOutputs/H1N1_A_Brisbane_59_2007_snp/H1N1_A_Brisbane_59_2007.fasta' #alex comment:same as bwa.
#opt$refGenomeGenBank <- '/home/agmcfarland/flu_project/FluPipeline/run_test/data/H1N1_A_Brisbane_59_2007.gb'
#opt$bcftoolsBin <- '/home/agmcfarland/miniconda3/envs/testenv/bin'
#opt$samtoolsBin <- '/home/agmcfarland/miniconda3/envs/testenv/bin'
#opt$bwaPath <- '/home/agmcfarland/miniconda3/envs/testenv/bin/bwa'

#opt$softwareDir <- '/home/agmcfarland/flu_project/FluPipeline'
#opt$process_method <- 'bushman_artic_v2'
#opt$outputFile <- '/home/agmcfarland/flu_project/shared_data/test_1_sample2/sampleOutputs/ashley_5/ashley_5.Rdata'
#opt$workDir <- '/home/agmcfarland/flu_project/shared_data/test_1_sample2/sampleOutputs/ashley_5'
#opt$R1 <- '/home/agmcfarland/flu_project/shared_data/test_data_1_sample2/ashley_5_R1_001.fastq.gz'
#opt$R2 <- '/home/agmcfarland/flu_project/shared_data/test_data_1_sample2/ashley_5_R2_001.fastq.gz'
#opt$refGenomeBWA <- '/home/agmcfarland/flu_project/shared_data/test_1_sample2/sampleOutputs/ashley_5/H1N1pdm_ref.fasta'
#opt$refGenomeFasta <- '/home/agmcfarland/flu_project/shared_data/test_1_sample2/sampleOutputs/ashley_5/H1N1pdm_ref.fasta' #alex comment:same as bwa.
#opt$refGenomeGenBank <- '/home/agmcfarland/flu_project/FluPipeline/references/H1N1pdm_ref.gb'
#opt$bcftoolsBin <- '/home/agmcfarland/miniconda3/envs/testenv/bin'
#opt$samtoolsBin <- '/home/agmcfarland/miniconda3/envs/testenv/bin'
#opt$bwaPath <- '/home/agmcfarland/miniconda3/envs/testenv/bin/bwa'

#opt$softwareDir <- '/home/agmcfarland/flu_project/FluPipeline'
#opt$process_method <- 'bushman_artic_v2'
#opt$outputFile <- '/home/agmcfarland/test_flu/output/sampleOutputs/HUP_19_007_S3/HUP_19_007_S3.Rdata'
#opt$workDir <- '/home/agmcfarland/test_flu/output/sampleOutputs/HUP_19_007_S3'
#opt$R1 <- '/home/agmcfarland/test_flu/data/HUP_19_007_S3_R1_001.fastq.gz'
#opt$R2 <- '/home/agmcfarland/test_flu/data/HUP_19_007_S3_R2_001.fastq.gz'
#opt$refGenomeBWA <- '/home/agmcfarland/test_flu/output/sampleOutputs/HUP_19_007_S3/H1N1pdm_ref.fasta'
#opt$refGenomeFasta <- '/home/agmcfarland/test_flu/output/sampleOutputs/HUP_19_007_S3/H1N1pdm_ref.fasta' #alex comment:same as bwa.
#opt$refGenomeGenBank <- '/home/agmcfarland/flu_project/FluPipeline/references/H1N1pdm_ref.gb'
#opt$bcftoolsBin <- '/home/agmcfarland/miniconda3/envs/testenv/bin'
#opt$samtoolsBin <- '/home/agmcfarland/miniconda3/envs/testenv/bin'
#opt$bwaPath <- '/home/agmcfarland/miniconda3/envs/testenv/bin/bwa'
#alex end:FOR TESTING MANUALLY assign opt vals

#alex start:set working directory to workDir
setwd(opt$workDir)
#alex end:set working directory to workDir

#alex start:copy bwa files to workDir and change the reference path locations to the files in workDir
system(paste0('cp ',opt$refGenomeBWA, ' .'))
system(paste0('cp ',opt$refGenomeFasta, ' .'))
opt$refGenomeBWA <- paste0(opt$workDir,'/',basename(opt$refGenomeBWA))
opt$refGenomeFasta <- paste0(opt$workDir,'/',basename(opt$refGenomeFasta))
#alex end:copy bwa files to workDir and change the reference path locations to the files in workDir

#alex start:MANUAL change source
source(paste0(opt$softwareDir, '/lib/assemble.lib.R'))#source(paste0(opt$softwareDir, '/lib/assemble.lib.R'))
#alex end:change source


if(opt$process_method == 'bushman_artic_v2'){
  opt$removeNTsFromAlignmentEnds <- 3
  opt$removeClippedReads <- FALSE #alex comment:changed this from TRUE to FALSE
}


# Handle missing required parameters.
if(! 'R1' %in% names(opt)) stop('--R1 must be defined.')
if(! 'R2' %in% names(opt)) stop('--R2 must be defined.')


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


#t1 <- paste0(opt$workDir, '/x') #alex comment: hashed out by alex. This makes the filename prefix `x`
#alex start: make the file prefix equal to the workdir's name.
samplename <- basename(opt$workDir)
t1 <- paste0(opt$workDir, '/', samplename)
#alex end:make the file prefix equal to the workdir's name.


# Combine R1 and R2 data files in composite R1 and R2 files.
system(paste0('cat ', paste0(R1s, collapse = ' '), ' > ', t1, '_R1.fastq'))
system(paste0('cat ', paste0(R2s, collapse = ' '), ' > ', t1, '_R2.fastq'))


# Quality trim reads and create trimmed FASTA files.
message('quality trimming...')
r <- prepareTrimmedReads(readFastq(paste0(t1, '_R1.fastq')), readFastq(paste0(t1, '_R2.fastq')), qualCode = opt$trimQualCode)
writeFasta(r[[1]], file = paste0(t1, '_R1.trimmed.fasta'))
writeFasta(r[[2]], file = paste0(t1, '_R2.trimmed.fasta'))


# Align trimmed reads to the reference genome.
#alex start:make bwa index file for reference
message('make bwa index...')
system(paste0(opt$bwaPath, ' index ', opt$refGenomeBWA)) #makes index files using the whole filename as the prefix
#alex end:make bwa index file for reference
message('align reads to reference...')
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



# Trim NTs from the ends of all aligned reads.#alex comment:this step is suggested by John to improve SNP detection, espcially with deletions
#alex comment: it doesn't seem like there is an option to not do this either.
#-------------------------------------------------------------------------------------------------- 
if(opt$removeNTsFromAlignmentEnds > 0){
  message('trim alignments...')
  message('Trimming alignments by ', opt$removeNTsFromAlignmentEnds > 0, ' NTs\n')
  system(paste0(opt$trimBamCommand, ' ', t1, '_genome.filt.bam ', t1, '_genome.filt.endTrim.bam ', opt$removeNTsFromAlignmentEnds, ' --clip'))
  system(paste0('cp ', t1, '_genome.filt.endTrim.bam ', t1, '_genome.filt.bam'))
}


# Remove read pairs with mapping qualities below the provided min score.
system(paste0(opt$samtoolsBin, '/samtools view -q ', opt$minBWAmappingScore, ' -b ', t1, '_genome.filt.bam > ', t1, '_genome.filt.qual.bam'))


# Retrieve a list of aligned reads.
alignedReadsIDs <- system(paste0(opt$samtoolsBin, '/samtools view ', t1, '_genome.filt.qual.bam | cut  -f 1 | uniq'), intern = TRUE)
opt$alignedReadIDs <- alignedReadsIDs

if(length(alignedReadsIDs) == 0){
  opt$errorCode <- 1
  opt$errorMessage <- 'No quality trimmed reads aligned to the reference genome.'
  save(opt, file = opt$outputFile)
  #unlink(opt$workDir, recursive = TRUE) #alex comment:hashed out by alex.
  stop()
}


# Sort and index filtered genome alignment.
message('sort and index genome alignment...')
system(paste0(opt$samtoolsBin, '/samtools sort -o ', paste0(t1, '_genome.filt.qual.sorted.bam'), ' ', paste0(t1, '_genome.filt.qual.bam')))
system(paste0(opt$samtoolsBin, '/samtools index ', paste0(t1, '_genome.filt.qual.sorted.bam')))


# Determine the maximum read depth. Overlapping mates will count a position twice.
#alex comment: changed the -o to > for the conda version of samtools. it is still v1.7 like John's version but it doesn't take -o for some reason.
#system(paste0(opt$samtoolsBin, '/samtools depth -d 0 ', paste0(t1, '_genome.filt.qual.sorted.bam'), ' -o ', paste0(t1, '.depth'))) #alex comment:depth file shows all segments have per base depth reported
message('calculate depth with samtools...')
system(paste0(opt$samtoolsBin, '/samtools depth -d 0 ', paste0(t1, '_genome.filt.qual.sorted.bam'), ' > ', paste0(t1, '.depth'))) #alex comment:depth file shows all segments have per base depth reported
maxReadDepth <- max(read.table(paste0(t1, '.depth'), sep = '\t', header = FALSE)[,3])

# Create pileup data file for determining depth at specific locations.
# (!) mpileup will remove duplicate reads.
#alex comment:http://www.htslib.org/doc/samtools-mpileup.html
#alex comment:-A(count orphans) -a(output all positions, even those with 0 depth) -Q(minimum base quality for a base to be considered) -d(max read depth to filter at)
message('pileup reads with samtools...')
system(paste0(opt$samtoolsBin, '/samtools mpileup -A -a -Q 0 -o ', paste0(t1, '.pileup'), ' -d ', maxReadDepth, 
              ' -f ', opt$refGenomeFasta, ' ', paste0(t1, '_genome.filt.qual.sorted.bam')))


# Determine the percentage of the reference genome covered in the pileup data.
#alex comment:don't try to print this variable in rstudio, it is huge
opt$pileupData <- tryCatch({
  read.table(paste0(t1, '.pileup'), sep = '\t', header = FALSE, quote = '')[,1:5]
}, error = function(e) {
  return(data.frame())
})


## Making variant table from pileup data ----------------------------------------------------------
# for testing
#setwd('/home/agmcfarland/flu_project/FluPipeline/run_test/output/sampleOutputs/H1N1_A_Brisbane_59_2007_snp')
#load('H1N1_A_Brisbane_59_2007_snp.Rdata')

# Pileup format reports the number of read pairs (column 4) while VCF format (DP) 
# reports the number of reads which appears to report 2x the pileup format value. 
# Confirmed by looking at pileup in IGV.

if(nrow(opt$pileupData) > 0){
  message('find all variants...')
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
}


##-----Find variants in pileup and vcf outputs. Adding to variant table #------------

if(nrow(opt$variantTable) > 0){
    message('Find variants in pileup and vcf outputs...')
    opt$variantTable <- tryCatch({
      # Here we parse the pileup data to create a more informative alt call for variants.
      #alex comment:for each row in the variant table generated by `bcftools call` and reference that row in the pileup data generated by `samtools pileup`.
      #alex comment:In the pileup data extract the string from the column containing the per read base stats http://www.htslib.org/doc/samtools-mpileup.html (keyword: This encodes information matches,)
      #alex comment:parsePileUpString2 counts the number of base substitutions observed uses grepl or indel info using a boolean flag to count indel length and sequence. It's a little complicated.
      x <- bind_rows(lapply(1:nrow(opt$variantTable), function(i){
        #browser()
      x <- opt$variantTable[i,]
      
      ###if(x$POS == 25445) browser()
      #d <- subset(opt$pileupData, V2 == x$POS) #original code hashed out by alex
      #alex comment: editing for multi chromosome support
      # subset the pileup data for the chromosome and position detected in the variant table
      d <- subset(opt$pileupData, V2 == x$POS & V1 == as.character(x$CHROM))
      #--------------------------------------
      # alex comment: search the d$V5 column which has indel, matching (.), and SNP [ACTG], as characters. The total number of each character equals its abundance at that position.
      # more info on pileup read data d$V5 column here: https://en.wikipedia.org/wiki/Pileup_format
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
    #browser()
    #opt$variantTable <- bind_rows(lapply(split(opt$variantTable, opt$variantTable$POS), function(x){ #original code hashed out by alex
    #browser() #original code hashed out by alex
    #x$reads <- subset(opt$pileupData, V2 == x$POS[1])$V4 #original code hashed out by alex
    #x$percentAlt <- x$altCount / x$reads[1] #original code hashed out by alex
    #x #original code hashed out by alex
    #})) #original code hashed out by alex
    #start new
    message('Merge variant table with pileup table...')
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
    
  } else {
    opt$variantTableMajor <- data.frame()
  }

###---------- Generate consensus sequence from filtered VCF and variantTableMajor #-------------------

# Generate consensus sequence from a filtered VCF file. Use the variantTable to filter for only major variants.
if(nrow(opt$variantTableMajor) > 0){
  message('Generate conensus sequence...')
  
  # Here we copy the variant vcf file and selectively remove calls 
  # which are not found in our filtered variantTable. This is done 
  # to preserve the additional comment lines which appear to be necessary. 
  
  tryCatch({
    # copy variant file and gunzip them 
    system(paste0('cp ', t1, '.filt.vcf.gz ', t1, '.filt.vcf.copy.gz'))
    system(paste0('gunzip ', t1, '.filt.vcf.copy.gz'))
    
    
    store_vcf_pass<- c()
    for (i in readLines(paste0(t1, '.filt.vcf.copy'))) {
      # for each line in filt.vcf.copy:
      # if the line starts with # then it is a header and will be stored
      if(grepl('^#', i)){
        store_vcf_pass <- c(store_vcf_pass,i)
        
      } else {
        # if the line does not start with #, then it is part of the vcf table. split it on \t and unlist it and store it in vector
        vcf_row <- unlist(strsplit(i,'\t'))
        
        # check if the vcf row is in the variantTableMajor
        df_puprow <- opt$variantTableMajor%>%filter(CHROM==vcf_row[1],POS==vcf_row[2])
        
        if (nrow(df_puprow) == 0){
          next
        } else {
          store_vcf_pass <- c(store_vcf_pass,i)
        }
      } 
    }
    
    write(unlist(store_vcf_pass[! sapply(store_vcf_pass, is.null)]), file = paste0(t1, '.filt.vcf.copy2'))
    
    # Capture vcf
    #read.table(textConnection(unlist(opt$finalVCF)), sep = '\t') #alex comment:hashed out by alex. unsure of what this does.
    #opt$finalVCF = o[! sapply(o, is.null)] #alex comment:hashed out by alex. unsure of what this does.
    
    # reformat the filtered vcf file into a gzip file
    system(paste0('bgzip ', t1, '.filt.vcf.copy2'))
    # index the filtered vcf file
    system(paste0(opt$bcftoolsBin, '/bcftools index ', t1, '.filt.vcf.copy2.gz'))
    # extract the consensus sequence using bcftools
    system(paste0('cat  ', opt$refGenomeFasta, ' | ', opt$bcftoolsBin, '/bcftools consensus ',  
                  t1, '.filt.vcf.copy2.gz > ', t1, '.consensus.fasta'))
    
    opt$concensusSeq <- readFasta(paste0(t1, '.consensus.fasta'))
    
  }, error=function(cond) {
    return('Error creating concensus sequence from vcf ')
  })
  
# if the major variant table does not have any variants then the consensus sequence is equivalent to the reference sequence. 
} else {
  opt$concensusSeq <- readFasta(opt$refGenomeFasta)
}

# write the fasta to file
writeFasta(opt$concensusSeq, paste0(t1, '.consensus.fasta')) 


# BCFtools calls indels by the base preceding the modification.
# Here we find deletion calls and increment the position by one to mark the first deleted base.
i <- opt$variantTable$POS %in% opt$variantTable[grep('del', opt$variantTable$ALT),]$POS
if(any(i)) opt$variantTable[i,]$POS <- opt$variantTable[i,]$POS+1

i <- opt$variantTableMajor$POS %in% opt$variantTableMajor[grep('del', opt$variantTableMajor$ALT),]$POS
if(any(i)) opt$variantTableMajor[i,]$POS <- opt$variantTableMajor[i,]$POS+1

message('write variant tables to file...')
write.csv(opt$variantTable,'variantTable.csv',row.names=FALSE)
write.csv(opt$variantTableMajor,'variantTableMajor.csv',row.names=FALSE)

message('write final rdata object to file...')
save(opt, file = opt$outputFile)
message('end of assemble.R')

