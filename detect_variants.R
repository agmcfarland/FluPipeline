
library(dplyr)
library(stringr)
library(optparse)
library(tidyr)

message('Beginning variant calling and consensus sequence generation...')

# # All options default to the default process_method.
option_list = list(
  make_option(c("--softwareDir"), type="character", default=NULL, help="path to R file", metavar="character"),
  make_option(c("--workDir"), type="character", default=NULL, help="path to work directory", metavar="character"),
  make_option(c("--R1"), type="character", default=NULL, help="comma delimited list of R1 fastq files", metavar="character"),
  make_option(c("--R2"), type="character", default=NULL, help="comma delimited list of R2 fastq files", metavar="character"),
  make_option(c("--refGenomeFasta"), type="character", default=NULL, help="reference genome FASTA path", metavar="character"),
  make_option(c("--minAmpliconLength"), type="integer", default=50, help="comma delimited list of R2 fastq files", metavar="character"),
  make_option(c("--maxAmpliconLength"), type="integer", default=350, help="comma delimited list of R2 fastq files", metavar="character"),
  make_option(c("--minVariantPhredScore"), type="integer", default=12, help="minimum PHRED score allowed for called varinats", metavar="character"),
  make_option(c("--removeNTsFromAlignmentEnds"), type="integer", default=3,  help="Number of NTs to remove from the end of alignments", metavar="character"),
  make_option(c("--BWAmappingScore"), type="integer", default=30, help="minimum BWA mapping score", metavar="character")
  )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

detect_VariantsMessages <- c('No reads were aligned','No pileup was created','No variants were detected')
detect_variants_message_file <- 'detect_variants_messages.txt'

## alex start:FOR TESTING MANUALLY assign opt vals
## use /home/agmcfarland/test_flu/output to store outputs
# opt$softwareDir <- '/home/agmcfarland/flu_project/FluPipeline'
# opt$workDir <- '/home/agmcfarland/quick_tests/output/sampleOutputs/Ashley_1_2_S81' #'/home/agmcfarland/flu_project/FluPipeline/run_test/output/sampleOutputs/IBV_Yamagata_ref_snp'
# opt$R1 <- 'fastp_trimmed_Ashley_1_2_S81_R1_001.fastq.gz' #'fastp_trimmed_IBV_Yamagata_ref_snp_R1_001.fastq.gz'
# opt$R2 <- 'fastp_trimmed_Ashley_1_2_S81_R2_001.fastq.gz' #'fastp_trimmed_IBV_Yamagata_ref_snp_R2_001.fastq.gz'
# opt$refGenomeFasta <- '/home/agmcfarland/flu_project/FluPipeline/references/H3N2_ref.fasta' #'/home/agmcfarland/flu_project/FluPipeline/references/IBV_Yamagata_ref.fasta'

# opt$softwareDir <- '/home/agmcfarland/flu_project/FluPipeline'
# opt$workDir <- '/home/agmcfarland/quick_tests/output/sampleOutputs/H3N2_test'
# opt$R1 <- 'fastp_trimmed_H3N2_test_R1_001.fastq.gz'
# opt$R2 <- 'fastp_trimmed_H3N2_test_R2_001.fastq.gz'
# opt$refGenomeFasta <- '/home/agmcfarland/flu_project/FluPipeline/references/H3N2_ref.fasta' #'/home/agmcfarland/flu_project/FluPipeline/references/IBV_Yamagata_ref.fasta'

# opt$softwareDir <- '/home/agmcfarland/flu_project/FluPipeline'
# opt$workDir <- '/home/agmcfarland/flu_project/test_variant_genomes/output/sampleOutputs/H3N2_test/H3N2_test_ivar'
# opt$R1 <- '/home/agmcfarland/flu_project/test_variant_genomes/output/sampleOutputs/H3N2_test/fastp_trimmed_H3N2_test_R1_001.fastq.gz'
# opt$R2 <- '/home/agmcfarland/flu_project/test_variant_genomes/output/sampleOutputs/H3N2_test/fastp_trimmed_H3N2_test_R2_001.fastq.gz'
# opt$refGenomeFasta <- '/home/agmcfarland/flu_project/test_variant_genomes/output/sampleOutputs/H3N2_test/H3N2_test.consensus.fasta' #'/home/agmcfarland/flu_project/FluPipeline/references/IBV_Yamagata_ref.fasta'
# ## End

## End
# 

# Sys.setenv(PATH="/home/agmcfarland/miniconda3/envs/testenv2/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games")

# set working directory and load helper functions
setwd(opt$workDir)

# define sample name and specify explicit path
samplename <- basename(opt$workDir)

message('Making bwa index...')
# copy reference genome fasta to the working directory to use for indexing
system(paste0('cp ',opt$refGenomeFasta, ' .'))

# get basename so indexing is done in the working directory
opt$refGenomeFasta <- basename(opt$refGenomeFasta)

# make index
system(paste0('bwa index ', opt$refGenomeFasta)) #makes index files using the whole filename as the prefix

# align reads to reference. make a sam file
message('Aligning reads to reference...')
system(paste0('bwa mem -M ', opt$refGenomeFasta, ' ',  opt$R1, ' ', opt$R2,  ' > ', paste0(samplename, '_genome.sam')))
gc()

# convert sam to bam and delete sam file (save space!)
system(paste0('samtools view -S -b -h ', paste0(samplename, '_genome.sam'), ' > ', paste0(samplename, '_genome.bam')))
invisible(file.remove(paste0(samplename, '_genome.sam')))

# Remove non-paired mates and filter alignment by insert size
system(paste0('samtools view -F 8 -h ',paste0(samplename,'_genome.bam'), " | awk 'length($10) > ", opt$minAmpliconLength, " && length($10) < ", opt$maxAmpliconLength, " || $1 ~ /^@/' | ", 'samtools view -b -h > ',paste0(samplename,'_genome.filt.bam')))

# trim ends of alignments
# use trimBam command and then rename the output file to be suffixed with _genome.filt.bam
if(opt$removeNTsFromAlignmentEnds > 0){
  message('Trimming alignments...')
  system(paste0('bam trimBam ', samplename, '_genome.filt.bam ', samplename, '_genome.filt.endTrim.bam ', opt$removeNTsFromAlignmentEnds, ' --clip'))
  system(paste0('cp ', samplename, '_genome.filt.endTrim.bam ', samplename, '_genome.filt.bam')) # copy over the original .filt.bam file
  invisible(file.remove(paste0(samplename, '_genome.filt.endTrim.bam'))) #remove to cleanup
}

# Remove read pairs with mapping qualities below the provided min score.
system(paste0('samtools view -h -q ', opt$BWAmappingScore, ' -b ', samplename, '_genome.filt.bam > ', samplename, '_genome.filt.qual.bam'))

# Get number of aligned reads and check that it is greater than 0. if it is not then write a file with the error 
lengthAlignedReadsIDs <- length(system(paste0('samtools view ', samplename, '_genome.filt.qual.bam | cut -f 1 | uniq'), intern = TRUE))

# exit the program if no aligned reads were found.
if(lengthAlignedReadsIDs == 0){
  message(detect_VariantsMessages[1])
  fileConn<-file(detect_variants_message_file)
  writeLines(c(detect_VariantsMessages[1]), fileConn)
  close(fileConn)
  stop()
}

message(paste0(lengthAlignedReadsIDs),' reads aligned')

# Sort and index filtered genome alignment.
message('Sorting and and indexing genome alignment...')
system(paste0('samtools sort -o ', paste0(samplename, '_genome.filt.qual.sorted.bam'), ' ', paste0(samplename, '_genome.filt.qual.bam')))
system(paste0('samtools index ', paste0(samplename, '_genome.filt.qual.sorted.bam')))

# Pileup reads
message('Piling up reads to get coverage...')
system(paste0('pileup.sh in=',paste0(samplename, '_genome.filt.qual.sorted.bam '), 'ref=',opt$refGenomeFasta,' out=',paste0(samplename,'_covstats.txt '),'basecov=',paste0(samplename,'_basecov.txt '),'countgc=f overwrite=t'))

if (nrow(read.csv(paste0(samplename,'_covstats.txt'),sep='\t')) == 0 ) {
  message(detect_VariantsMessages[2])
  fileConn<-file(detect_variants_message_file)
  writeLines(c(detect_VariantsMessages[2]), fileConn)
  close(fileConn)
  stop()
}

message('Calling variants...')

system(paste0('callvariants.sh in=',paste0(samplename, '_genome.filt.qual.sorted.bam '), 'ref=',opt$refGenomeFasta,' out=',paste0(samplename, '_allVariants.vcf '), 'shist=',paste0(samplename, '_variantQualityHisto.txt '), 
              'rarity=0 overwrite=t clearfilters'))


# read in variant .vcf file
allVariantsRaw <- tryCatch({
  allVariantsRaw <- read.csv(paste0(samplename, '_allVariants.vcf'), comment.char = '#',sep='\t',header=FALSE)
}, error = function(e) {
  
  message(detect_VariantsMessages[3])
  fileConn<-file(detect_variants_message_file)
  writeLines(c(detect_VariantsMessages[3]), fileConn)
  close(fileConn)
  
  # write all empty dataframes for different variantTable levels
  write.csv(data.frame(),paste0(samplename,'_variantTableRaw.csv'),row.names=FALSE)
  write.csv(data.frame(),paste0(samplename,'_variantTable.csv'),row.names=FALSE)
  write.csv(data.frame(),paste0(samplename,'_variantTableMajor.csv'),row.names=FALSE)
  stop()
})
colnames(allVariantsRaw) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','GENOME')

# split the INFO column into several columns using the `;` character as the marker
allVariantsRaw$INFO <- as.character(allVariantsRaw$INFO)
split_colnames <- unlist(str_split(allVariantsRaw$INFO[1], ';'))
split_colnames <- str_split(split_colnames, '=')
split_colnames_final <- c()
for (i in split_colnames){
  split_colnames_final <- c(split_colnames_final,i[1])
}
allVariantsRaw <- allVariantsRaw%>%
  separate(INFO,sep=';',into=split_colnames_final)%>%
  mutate(across(everything(),~ gsub('.*=','', .)))

# calculate variant depth
allVariantsFiltered <- merge(read.csv(paste0(samplename,'_basecov.txt'),sep='\t')%>%rename('CHROM'=1,'POS'=2,'pileup_cov'=3),
                        allVariantsRaw, 
                        by = c('CHROM','POS'),all.y=TRUE)

allVariantsFiltered$QUAL <- as.numeric(allVariantsFiltered$QUAL)
allVariantsFiltered$AF <- as.numeric(allVariantsFiltered$AF)

# Create a filtered variant file
# opt$minVariantPhredScore <- 5
allVariantsFiltered <- allVariantsFiltered%>%filter(QUAL>=opt$minVariantPhredScore)

# get major variants
majorVariants <- allVariantsFiltered%>%
  filter(AF>=0.5)

message('write variant tables to file...')
write.csv(allVariantsRaw,paste0(samplename,'_variantTableRaw.csv'),row.names=FALSE)
write.csv(allVariantsFiltered,paste0(samplename,'_variantTable.csv'),row.names=FALSE)
write.csv(majorVariants,paste0(samplename,'_variantTableMajor.csv'),row.names=FALSE)

message('Finished variant calling')




