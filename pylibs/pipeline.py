import sys
import os 
from os.path import join as pjoin
import multiprocessing as mp
from multiprocessing import Pool
import pandas as pd
import re
import glob
import shutil
import subprocess
from datetime import datetime
from Bio.Seq import Seq
from Bio import SeqIO, SeqFeature
import argparse
import logging
import json
import numpy as np
from .processing_classes import SequencingSample, RunLogger
from .processing_functions import *
from .helper_functions import call_Command


def flu_Pipeline(
	baseDirectory,
	softwareDir,
	sample,
	sequenceDataDir,
	referenceStrainsDir,
	force_overwrite,
	keep_all_intermediate_files,
	strain_sample_depth,
	downsample,
	consensus_masking_threshold,
	min_variant_phred_score,
	remove_NTs_from_alignment_ends,
	min_read_mapping_score,
	masked_nextclade,
	masked_ivar,
	base_quality,
	keep_duplicates,
	min_variant_frequency,
	use_strain,
	keep_trimmed_reads,
	major_variant_frequency,
	major_indel_frequency,
	minimum_read_depth,
	major_variant_caller,
	intrahost_variant_caller,
	single_pass
	):
	'''
	Runs through all steps of the flu pipeline. 
	'''

# # ## troubleshooting inputs start ##
# ## use this exact input for troubleshooting the pipeline
# os.chdir('/data/flu_project/benchmarking_project/output')
# Rscript = 'Rscript'
# softwareDir = '/data/FluPipeline'
# script_path = '/data/FluPipeline'
# baseDirectory = '/data/flu_project/benchmarking_project/output'
# # referenceStrainsDir = '/data/FluPipeline/references'
# # sequenceDataDir = '/data/flu_project/benchmarking_project/data'
# referenceStrainsDir = '/data/flu_project/benchmarking_project/data/mccrone/references'
# sequenceDataDir = '/data/flu_project/benchmarking_project/data/mccrone/reads'
# use_fasta = False
# force_overwrite = True
# force_base_directory = True
# keep_all_intermediate_files = True
# threads = 8
# max_mem_per_thread = 4
# strain_sample_depth = 2000
# downsample = -1
# base_quality = 30
# no_deduplicate = False
# remove_NTs_from_alignment_ends = 3
# assembly = False
# min_read_mapping_score = 30
# min_variant_phred_score = 20
# min_variant_frequency = 0.05
# consensus_masking_threshold = 0
# masked_nextclade = False
# masked_ivar = False
# runtest = False
# variant_caller = 'bbtools'
# detect_strain = True
# use_strain = None

# sample = SequencingSample()
# sample.get_DataFromReadPairs(read1_filename=pjoin(sequenceDataDir,'SRR6121517_1_S1_R1_001.fastq.gz'))
# sample.get_SampleDirPath(directory=baseDirectory)
# sample.check_ReadExistence(directory=sequenceDataDir)
## troubleshooting inputs end ##

	## ==================================Prepare sample run directory================================================
	## ==============================================================================================================
	# set and create directory to store all intermediate data in sample-specific folders. Switch to that directory afterwards.
	run_chunk = True

	if force_overwrite==True:
		try: 
			shutil.rmtree(sample.dirpath)
		except:
			pass

	os.makedirs(sample.dirpath, exist_ok=True)
	os.chdir(sample.dirpath)
	intrahost_var_directory = pjoin(sample.dirpath,sample.samplename+'_ivar')
	os.makedirs(intrahost_var_directory,exist_ok=True)

	start_run_timer = datetime.now()
	# record bulk ouputs to sample_logger
	sample_logger = RunLogger(directory=pjoin(baseDirectory,'sampleLogs'),filename=sample.samplename)
	sample_logger.initialize_FileHandler()
	sample_logger.logger.info('Processing sample: {}'.format(sample.samplename))
	# use print so that the program talks more
	print('Processing sample: {}'.format(sample.samplename))

	## ==================================Find and assign reference strain============================================
	## ==============================================================================================================
	if run_chunk == True:
		try:
			sample_logger.logger.info('Assigning reference strain')

			# if no_deduplicate == False:
			# 	# read deduplication. quality, adaptor trimming with fastp
			# 	sample_logger.logger.info('Deduplicate reads, read trimming and quality scoring\n')
			# 	call_Command(cmd=
			# 		[
			# 		'fastp', 
			# 		'--in1', pjoin(sequenceDataDir,sample.read1_filename), #readfileR1
			# 		'--in2', pjoin(sequenceDataDir,sample.read2_filename), #readfileR2
			# 		'--out1', 'fastp_trimmed_{}'.format(sample.read1_filename), 
			# 		'--out2', 'fastp_trimmed_{}'.format(sample.read2_filename),
			# 		'-j', 'fastp_stats_{}.json'.format(sample.samplename), #json output
			# 		'-h', 'fastp_stats_{}.html'.format(sample.samplename), #json output
			# 		'-q', str(base_quality),#, #quality score 30
			# 		'overwrite=True',
			# 		'--cut_front', #cut from front
			# 		'--cut_window_size', '4',#, #cut_window_size
			# 		'--cut_mean_quality', '20',
			# 		'--length_required', '50', #minimum read length
			# 		'--thread', '3',
			# 		'--dedup'
			# 		],
			# 		logger_=sample_logger)

			# else :
			# quality, adaptor trimming with fastp
			sample_logger.logger.info('Read trimming and quality scoring\n')
			call_Command(cmd=
				[
				'fastp', 
				'--in1', pjoin(sequenceDataDir,sample
				.read1_filename), #readfileR1
				'--in2', pjoin(sequenceDataDir,sample.read2_filename), #readfileR2
				'--out1', 'fastp_trimmed_{}'.format(sample.read1_filename), 
				'--out2', 'fastp_trimmed_{}'.format(sample.read2_filename),
				'-j', 'fastp_stats_{}.json'.format(sample.samplename), #json output
				'-h', 'fastp_stats_{}.html'.format(sample.samplename), #json output
				'-q', str(base_quality),#, #quality score 30
				'overwrite=True',
				'--cut_front', #cut from front
				'--cut_window_size', '4',#, #cut_window_size
				'--cut_mean_quality', '20',
				'--length_required', '50', #minimum read length
				'--thread', '3'
				],
				logger_=sample_logger)


			if downsample != -1:
				# downsamples to the user-inputted number of reads
				sample_logger.logger.info('Downsampling reads to {}\n'.format(downsample))
				call_Command(cmd=
					[		
					'reformat.sh', 
					'in1=fastp_trimmed_{}'.format(sample.read1_filename),
					'in2=fastp_trimmed_{}'.format(sample.read2_filename),
					'out1=fastp_trimmed1_{}'.format(sample.read1_filename),
					'out2=fastp_trimmed1_{}'.format(sample.read2_filename),
					'reads={}'.format(downsample),
					'overwrite=t'
					],
					logger_=sample_logger)
				os.system('cp fastp_trimmed1_{} fastp_trimmed_{}'.format(sample.read1_filename, sample.read1_filename))
				os.system('cp fastp_trimmed1_{} fastp_trimmed_{}'.format(sample.read2_filename, sample.read2_filename))
				os.remove('fastp_trimmed1_{}'.format(sample.read1_filename))
				os.remove('fastp_trimmed1_{}'.format(sample.read2_filename))


			# convert fastp json into table
			with open('fastp_stats_{}.json'.format(sample.samplename)) as infile:
				readstats = json.load(infile)
			pd.DataFrame(
				{
				'sample':[sample.samplename],
				'total_reads_before':[readstats['summary']['before_filtering']['total_reads']],
				'total_reads_after':[readstats['summary']['after_filtering']['total_reads']],
				'passed_filter_reads':[readstats['filtering_result']['passed_filter_reads']],
				'low_quality_reads':[readstats['filtering_result']['low_quality_reads']],
				'too_many_N_reads':[readstats['filtering_result']['too_many_N_reads']]
				}).to_csv('fastp_stats_{}.csv'.format(sample.samplename), index=None)

			if use_strain == None:
				# randomly select a subset of reads with reformat.sh
				sample_logger.logger.info('Subsample reads\n')
				call_Command(cmd=
					[
					'reformat.sh',
					'in1=fastp_trimmed_{}'.format(sample.read1_filename), 
					'in2=fastp_trimmed_{}'.format(sample.read2_filename),
					'out1=subset_fastp_trimmed_{}'.format(sample.read1_filename), 
					'out2=subset_fastp_trimmed_{}'.format(sample.read2_filename),
					'samplereadstarget={}'.format(strain_sample_depth),
					'overwrite=True'
					],
					logger_=sample_logger)

				# find best reference
				os.chdir(softwareDir)
				sample_logger.logger.info('Summarize coverage and pick reference\n')
				call_Command(cmd=
					[
					'python', '-m', 'bin.best_reference',
					'--baseDir', '{}'.format(pjoin(sample.dirpath)),
					'--logDir', '{}'.format(pjoin(baseDirectory,'sampleLogs')),
					'--referenceStrainsDir', '{}'.format(referenceStrainsDir),
					'--R1', '{}/subset_fastp_trimmed_{}'.format(sample.dirpath,sample.read1_filename),
					'--R2', '{}/subset_fastp_trimmed_{}'.format(sample.dirpath,sample.read1_filename)
					],logger_=sample_logger
					)
				os.chdir(sample.dirpath)

				os.remove(f'subset_fastp_trimmed_{sample.read1_filename}')
				os.remove(f'subset_fastp_trimmed_{sample.read2_filename}')

				# get the best reference strain
				reference_strain = select_BestReference(samplename=sample.samplename, directory=sample.dirpath)
			else:
				reference_strain = use_strain
				if reference_strain.endswith('.gb') == True:
					reference_strain.replace('.gb','.fasta')
				if reference_strain.endswith('.fasta') == False:
					reference_strain = reference_strain+'.fasta'

		except:
			sample_logger.logger.exception('FluPipeline Error: failure at strain assignment', exc_info=True)
			run_chunk = False

	
	# set new trimmed filenames
	trimmed_r1 = 'fastp_trimmed_'+sample.read1_filename
	trimmed_r2 = 'fastp_trimmed_'+sample.read2_filename
	
	# set duplicate read flag
	if keep_duplicates == True:
		keep_duplicate_flag = '--keep_duplicates'
	else:
		keep_duplicate_flag = ''

	## ==================================find variants===============================================================
	## ==============================================================================================================
	if run_chunk == True:
		try:
			sample_logger.logger.info('Finding variants')

			os.chdir(softwareDir)
			call_Command(cmd=
				[
				'python', '-m', 'bin.detect_variants',
				'--baseDir', '{}'.format(pjoin(sample.dirpath)),
				'--logDir', '{}'.format(pjoin(baseDirectory,'sampleLogs')),
				'--variant_caller', '{}'.format(pjoin(major_variant_caller)),
				'--R1', '{}'.format(pjoin(sample.dirpath,trimmed_r1)),
				'--R2', '{}'.format(pjoin(sample.dirpath,trimmed_r2)),
				'--refGenomeFasta', '{}'.format(pjoin(referenceStrainsDir,reference_strain)),
				'--minVariantPhredScore', '{}'.format(min_variant_phred_score),
				'--removeNTsFromAlignmentEnds', '{}'.format(remove_NTs_from_alignment_ends),
				'--BWAmappingScore', '{}'.format(min_read_mapping_score),
				'--minorVariantThreshold', '{}'.format(min_variant_frequency),
				'--consensus_sequence',
				'--consensus_masking_threshold', '{}'.format(consensus_masking_threshold),
				'--majorIndelVariantThreshold', '{}'.format(major_variant_frequency),
				'--majorVariantThreshold', '{}'.format(major_indel_frequency),
				'--minimum_read_depth', '{}'.format(minimum_read_depth)
				],logger_=sample_logger
				)
			os.chdir(sample.dirpath)

		except:
			sample_logger.logger.exception('FluPipeline Error: failure at variant calling', exc_info=True)
			run_chunk = False

	## ==================================nextclade clade assignment==================================================
	## ==============================================================================================================
	if run_chunk == True:
		try:
			if masked_nextclade == True:
				# nextclade_input_fasta = sample.samplename+'_masked_consensus_sequence.fasta'
				nextclade_input_fasta = sample.samplename+'_consensus_sequence.fasta'
			else:
				nextclade_input_fasta = sample.samplename+'_consensus_sequence.fasta'

			# check whether the nextclade reference is available for the given reference strain used by the sample
			nextclade_reference = reference_NextCladeLookUp(reference=reference_strain.replace('.fasta',''))

			if nextclade_reference != 'none available':
				sample_logger.logger.info('Finding clade assignment using nextclade')

				# make new fasta file containing only the HA sequence (segment 4) that is used by nextclade for clade assignment
				with open('{}_HA.consensus.fasta'.format(sample.samplename),'w') as infile:
					for record in SeqIO.parse(nextclade_input_fasta,'fasta'):
						if record.id.endswith('_4') == True:
							infile.write('>{}-{}\n'.format(sample.samplename,record.id))
							infile.write('{}\n'.format(str(record.seq)))

				# copy the reference directory files from the nextclade_references directory and store it in the sample directory			
				shutil.copytree(src=pjoin(softwareDir,'nextclade_references',nextclade_reference), dst=pjoin(os.getcwd(),nextclade_reference))

				# run nextclade. all ouputs are stored in a folder called 'nextclade_output' created by the nextclade program
				call_Command(cmd=
					[
					'nextclade',
					'run', 
					'--in-order',
					'--output-basename', sample.samplename,
					'--input-dataset', pjoin(os.getcwd(),nextclade_reference),
					'--output-all=nextclade_output/',
					'{}_HA.consensus.fasta'.format(sample.samplename)
					],
					logger_=sample_logger)

				shutil.rmtree(nextclade_reference)

			else:
				sample_logger.logger.info('No nextclade reference strain exists for this sample strain')

		except:
			sample_logger.logger.exception('FluPipeline Error: failure at nextclade clade assignment', exc_info=True)
			run_chunk = False

	
	## ==================================Intra-host variant detection================================================
	## ==============================================================================================================
	if run_chunk == True:
		if single_pass == False:
			try:
				sample_logger.logger.info('Calling variants from consensus sequence')

				if masked_ivar == True:
					# ivar_input_fasta = sample.samplename+'_masked_consensus_sequence.fasta'
					ivar_input_fasta = sample.samplename+'_consensus_sequence.fasta'
				else:
					ivar_input_fasta = sample.samplename+'_consensus_sequence.fasta'

				os.chdir(softwareDir)
				call_Command(cmd=
					[
					'python', '-m', 'bin.detect_variants',
					'--baseDir', '{}'.format(pjoin(intrahost_var_directory)),
					'--logDir', '{}'.format(pjoin(baseDirectory,'sampleLogs')),
					'--variant_caller', '{}'.format(pjoin(intrahost_variant_caller)),
					'--R1', '{}'.format(pjoin(sample.dirpath,trimmed_r1)),
					'--R2', '{}'.format(pjoin(sample.dirpath,trimmed_r2)),
					'--refGenomeFasta', '{}'.format(pjoin(sample.dirpath,ivar_input_fasta)),
					'--minVariantPhredScore', '{}'.format(min_variant_phred_score),
					'--removeNTsFromAlignmentEnds', '{}'.format(remove_NTs_from_alignment_ends),
					'--BWAmappingScore', '{}'.format(min_read_mapping_score),
					'--minorVariantThreshold', '{}'.format(min_variant_frequency),
					'--consensus_masking_threshold', '{}'.format(consensus_masking_threshold),
					'--majorIndelVariantThreshold', '{}'.format(major_variant_frequency),
					'--majorVariantThreshold', '{}'.format(major_indel_frequency),
					'--minimum_read_depth', '{}'.format(minimum_read_depth)
					],logger_=sample_logger
					)
				os.chdir(sample.dirpath)

			except:
				sample_logger.logger.exception('FluPipeline Error: failure at intrahost variant calling', exc_info=True)
				sample_logger.logger.info('FluPipeline Time: {}'.format(datetime.now()-start_run_timer))
				run_chunk = False

	# ==================================clean up and record run=====================================================
	# ==============================================================================================================
	if keep_trimmed_reads == False:
		os.remove(trimmed_r1)
		os.remove(trimmed_r2)

	if run_chunk == True:
		sample_logger.logger.info('Finished processing sample {}'.format(sample.samplename))
	else:
		sample_logger.logger.info('Finished processing sample {} with warnings'.format(sample.samplename))

	sample_logger.logger.info('FluPipeline Time: {} (H:M:S)'.format(datetime.now()-start_run_timer))
	#use print so the program talks more
	print('Finished processing sample {}'.format(sample.samplename))
	print('FluPipeline Time: {}'.format(datetime.now()-start_run_timer))

