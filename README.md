# FluPipeline

`FluPipeline` is a command line program for processing influenza sequencing data. It takes as input a folder of paired-end read files and outputs a folder containing variant call, taxonomic, and mapping data and visualizations.

## Overview

<img src="docs/overview.png" alt="FluPipeline overview" width=1000>

**Highlighted Features**

- Detailed tables and figures
- Custom reference strains can be added
- Automatic thread selection selection to run samples in parallel
- Downsample reads and run pipeline
- Automatic testing
- Simple installation

---

## Installation instructions

Download the FluPipeline repository

Follow the directions below to create a conda environment with all required packages and dependencies.

```
cd /path/to/FluPipeline-main/install

yes | conda env create --name FluPipeline_env --file environment.yml

conda activate FluPipeline_env

yes | Rscript ./r_packages.R
```

## Run a test

Generates synthetic read data and runs FluPipeline on it. Great for making sure everything installed correctly!

```
cd /path/to/FluPipeline-main

conda activate FluPipeline_env

python FluPipeline.py \
--runtest \
--threads 6
```

## Usage examples

Run FluPipeline on fastq files in a folder (--sequence_directory) and outputs to a new folder (--base_directory). Removes the output folder if it already exists (--force_base_directory). Chooses the appropritate number of threads to run samples in parallel based on the assumption one sample will use 8 Gb of memory (--max_mem_per_thread).

### Basic

```
python FluPipeline.py \
--base_directory /path/to/output/folder \
--sequence_directory /path/to/sequence/folder \
--max_mem_per_thread 8 \
--force_base_directory
```

The same as above but without the de-novo assembly step (--no_assembly).

```
python FluPipeline.py \
--base_directory /path/to/output/folder \
--sequence_directory /path/to/sequence/folder \
--max_mem_per_thread 8 \
--force_base_directory \
--no_assembly
```


### Advanced

```
python FluPipeline.py \
--base_directory /path/to/output/folder \
--sequence_directory /path/to/sequence/folder \
--force_base_directory \ # overwrite base_directory if it exists
--threads 6 \ # use 6 threads
--use_fasta \ # use custom reference genomes (in fasta format)
--reference_directory /path/to/custom/references \ # directory with custom reference genomes
--downsample 1000000 \ # downsample all read pairs to 1 million reads
--strain_sample_depth 10000 # use 10,000 reads for reference strain identification (default is 2,000)
```

## Outputs

The output folder will have the following tree:

```
.
├── runLog.txt
├── runStats.csv
├── runSummary.pdf
├── sampleLogs
├── sampleOutputs
├── sampleReports
└── softwareVersions.txt
```

###  Output descriptions

**Folders**

1. **sampleOutputs:** Each sample read pair (from here on termed a sample) has a folder with the sample's name containing all files created during processing.

2. **sampleResults:** Each sample report is copied here for quick access.

3. **sampleLogs:** A record of how the run for each sample went.


**Files**

1. **RunSummary.pdf:** 

- Reference strains used

- Overall read quality and coverage for each sample

- Errors encountered in a readable format

2. **runLog.txt:** All messages during run time

3. **runStats.csv:** Samples processed, their runtime, and errors encountered

4. **softwareVersions.txt:** The names and versions for all packages/tools/programs/libraries used by FluPipeline


## Usage Parameters

All usage parameters supported by FluPipeline. Default values are in brackets.

```
usage: FluPipeline [-h] [--base_directory] [--reference_directory]
                   [--sequence_directory] [--use_fasta] [--force]
                   [--force_base_directory] [--keep_all_intermediate_files]
                   [--threads] [--max_mem_per_thread] [--strain_sample_depth]
                   [--downsample] [--base_quality] [--no_deduplicate]
                   [--remove_NTs_from_alignment_ends]
                   [--min_read_mapping_score] [--min_variant_phred_score]
                   [--min_variant_frequency] [--consensus_masking_threshold]
                   [--masked_nextclade] [--masked_ivar] [--runtest]

optional arguments:
  -h, --help            show this help message and exit

  # main arguments
  --base_directory      directory that run will output data to.
                        [./FluPipeline_output]
  --reference_directory
                        directory containing reference strain files (.gb or
                        .fasta (see --use_fasta flag))
                        [script_path/references]
  --sequence_directory
                        directory containing fastq sequence files (.gz format)
                        [None]
  --use_fasta           fasta format: fasta file(s) contain all eight segments
                        sequences. All segments must have a single name (only
                        letters, numbers, and underscores. At the end of the
                        name there should be an underscore followed by the
                        segment number. Example: an_example_name_1. [False]
  
	# remove files/folders
  --force               overwrite existing sample files. [False]
  --force_base_directory
                        overwrite existing directory. [False
  --keep_all_intermediate_files
                        remove intermediate files. [False]
  
  # parallel jobs                      
  --threads             number of samples to process in paralell. one sample
                        is one read pair [4]
  --max_mem_per_thread
                        automatically determines the number of threads to use
                        based on memory per thread supplied (in Gb) [None]
  
  # number of reads to use
  --strain_sample_depth
                        number of random reads to use to determine strain
                        assignment. [2000]
  --downsample          downsample all read files to these many reads. [-1 (no
                        downsampling)]
  
  # read processing                    
  --base_quality        keep reads that have at least an average of this
                        phred-scaled value. [30]
  --no_deduplicate      do not conduct read deduplication. [False]
  --remove_NTs_from_alignment_ends
                        remove this many bases from the left and right of each
                        read prior to mapping. [3]
  
  # read mapping
  --min_read_mapping_score
                        keep reads that mapped above or eequal to this phred-
                        scaled value. [3]
  
  # assembly
  --no_assembly         Do not assemble reads into a draft genome. [False]

  # variant calling 
  --min_variant_phred_score
                        keep all variants above or equal to this phred-scaled
                        value. [5]
  --min_variant_frequency
                        keep all variants with allele frequencies above or
                        equal this value. [0.05]
  
  # consensus sequence generation and usage
  --consensus_masking_threshold
                        replace any nucleotides in the consensus sequence with
                        N if their depth falls below this number. [0]
  --masked_nextclade    use the masked consensus sequence fasta file for
                        nextclade clade assignment. [False]
  --masked_ivar         use the masked consensus sequence fasta file as the
                        reference genome for intrahost variation detection.
                        [False]
  
  # run test
  --runtest             run an in silico test to make sure FluPipeline is
                        working correctly. [False]
  ```

# Citation

[FluPipeline - Alex McFarland - github.com/agmcfarland/FluPipeline](https://github.com/agmcfarland/FluPipeline)

[Follow me on twitter! @alexmcfarland_](https://twitter.com/alexmcfarland_)

