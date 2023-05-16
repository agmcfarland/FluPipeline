# FluPipeline

`FluPipeline` is a command line program for processing Influenza sequencing data. It takes as input a folder of paired-end read files, finds the best flu strain to find variants against, and outputs a variant calls, mapping and coverage data, and taxonomic placement.

## Overview

<img src="docs/overview.png" alt="FluPipeline overview" width=1000>

**Highlighted Features**

- Detailed VCF tables in CSV format for several different variant callers (bcftools, bbmap, lofreq, freebayes)
- Call intra host variants
- Automatic reference strain selection
- Custom reference strains can be added
- Automatic thread selection selection to run samples in parallel
- Downsample reads and run pipeline
- Simple installation
- Automated testing

---

## Installation

Enter `FluPipeline` directory and install/load environment.

```
cd /path/to/FluPipeline

conda env create --file /install/environment.yml

conda activate testenv
```

## Test installation

Verify that expected outputs are created by `FluPipeline`.

```
cd /path/to/FluPipeline

python -m unittest tests.integration.test_detect_variants -v

python -m unittest tests.integration.test_best_reference -v

python FluPipeline.py --runtest
```

## Example set up

Simple skeleton code to run`FluPipeline` on a set of reads.

```
cd /path/to/FluPipeline

python FluPipeline.py \
--base_directory /path/to/output \
--sequence_directory /path/to/reads \
--threads 8 \
--force_base_directory \
--intrahost_variant_caller bbtools
```

## detect_variants

Run a variant caller with different parameters on a specific sample. Useful if you want to try different settings/variant callers and compare while cutting down on computational load.

```
cd /path/to/FluPipeline

python -m bin.detect_variants \
--build_input_from /path/to/sampleOutputs/sample \
--output_name my_output_folder \
--variant_caller lofreq \
--BWAmappingScore 10 \
--minVariantPhredScore 5 \
--minimum_read_depth 1
```

`detect_variants.py` is the variant calling workhorse. You can supply reads, an output file, and a reference file. Check out all options by doing:

```
python -m bin.detect_variants -h
```

## Monitor run

In a different terminal, check how many samples have finished/are running/are errored out. Sometimes samples with poor sequencing depth will error out, but this will not stop FluPipeline from processing other samples.

```
cd /path/to/FluPipeline

python bin/check_finished.py /path/to/base/directory
```


# Options

FluPipeline has many options, check them out below:

```
cd /path/to/FluPipeline

python FluPipeline.py -h
```

```
usage: FluPipeline v0.7.0 [-h] [--base_directory] [--reference_directory] [--sequence_directory] [--use_fasta] [--force] [--force_base_directory] [--threads] [--max_mem_per_thread] [--strain_sample_depth] [--downsample] [--use_strain]
                          [--base_quality] [--keep_duplicates] [--remove_NTs_from_alignment_ends] [--keep_trimmed_reads] [--min_read_mapping_score] [--gap_open_penalty] [--gap_extension_penalty] [--major_variant_caller]
                          [--intrahost_variant_caller] [--single_pass] [--min_variant_phred_score] [--min_variant_frequency] [--major_variant_frequency] [--major_indel_frequency] [--minimum_read_depth] [--runtest]

options:
  -h, --help            show this help message and exit
  --base_directory      directory that run will output data to [./FluPipeline_output]
  --reference_directory
                        directory containing reference strain files (.gb or .fasta (see --use_fasta flag if you want to use fasta files)) [script_path/references]
  --sequence_directory
                        directory containing fastq sequence files (.gz format) [None]
  --use_fasta           fasta format: fasta file(s) contain all eight segments sequences. All segments must have a single name (only letters, numbers, and 3 underscores. At the end of the name there should be an underscore followed by the
                        segment number. Example: an_example_1. [False]
  --force               overwrite existing sample files. [False]
  --force_base_directory
                        overwrite existing directory. [False
  --threads             number of samples to process in paralell. one sample is one read pair [4]
  --max_mem_per_thread
                        automatically determines the number of threads to use based on memory per thread supplied (in Gb) [None]
  --strain_sample_depth
                        number of random reads to use to determine strain assignment. [2000]
  --downsample          downsample all read files to these many reads. [-1 (no downsampling)]
  --use_strain          name of strain. No file extensions in name. [None]
  --base_quality        keep reads that have at least an average of this phred-scaled value. [30]
  --keep_duplicates     do not conduct read deduplication. [False]
  --remove_NTs_from_alignment_ends
                        remove this many bases from the left and right of each read prior to mapping. [3]
  --keep_trimmed_reads  keep trimmed reads used for analysis. [False]
  --min_read_mapping_score
                        keep reads that mapped above or equal to this MAPQ value. [10]
  --gap_open_penalty    The gap open penalty used by BWA during alignment. [6]
  --gap_extension_penalty
                        The gap extension penalty used by BWA during alignment. [1]
  --major_variant_caller
                        variant caller to use (bcftools, bbmap, lofreq, freebayes). [bcftools]
  --intrahost_variant_caller
                        intra host variant caller to use (bcftools, bbmap, lofreq, freebayes). [lofreq]
  --single_pass         only call variants one time (no variants from consensus sequence will be called). [False]
  --min_variant_phred_score
                        keep all variants above or equal to this phred-scaled value. [20]
  --min_variant_frequency
                        keep all variants with allele frequencies above or equal this value. [0.05]
  --major_variant_frequency
                        keep all major variants with allele frequencies above or equal this value. [0.5]
  --major_indel_frequency
                        keep all major indels with allele frequencies above or equal this value. [0.8]
  --minimum_read_depth
                        Mask/ignore all bases and variants at or below this read depth. [10]
  --runtest             Run an in silico test to make sure FluPipeline is working correctly. [False]
  ```

# Citation

[FluPipeline - Alex McFarland - github.com/agmcfarland/FluPipeline](https://github.com/agmcfarland/FluPipeline)

[Follow me on twitter! @alexmcfarland_](https://twitter.com/alexmcfarland_)

