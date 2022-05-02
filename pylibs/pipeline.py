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
from .processing_classes import SequencingSample, RunLogger, ConsensusSeq
from .processing_functions import *
from .helper_functions import call_Command

def flu_Pipeline(
	baseDirectory,
	Rscript,
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
	no_deduplicate,
	min_variant_frequency,
	no_assembly):
	'''
	Runs through all steps of the flu pipeline. 
	'''

# ## troubleshooting inputs start ##
## use this exact input for troubleshooting the pipeline
# baseDirectory='/home/agmcfarland/quick_tests/output'
# Rscript='Rscript'
# softwareDir = '/home/agmcfarland/flu_project/FluPipeline'
# sample = SequencingSample()
# sample.get_DataFromReadPairs(read1_filename='Ashley_1_2_S81_R1_001.fastq.gz')
# sample.get_SampleDirPath(directory=baseDirectory)
# sequenceDataDir = '/home/agmcfarland/quick_tests/data'
# referenceStrainsDir = '/home/agmcfarland/flu_project/FluPipeline/references' #directory reference strains to find the best match for a given readset
# force_overwrite = True #deletes and remakes the sample folder in sampleOutputs
# keep_all_intermediate_files = False
# strain_sample_depth = 3000
# downsample = 500000
# consensus_masking_threshold = 0
# min_variant_phred_score = 4
# remove_NTs_from_alignment_ends = 3
# min_read_mapping_score = 30
# masked_nextclade= False
# masked_ivar = False
# base_quality = 30
# deduplicate = True


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
	detect_variants_message_file = 'detect_variants_messages.txt'

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
			# copy reference fastas to sample directory 
			reference_fasta_list = glob.glob(pjoin(referenceStrainsDir,'*.fasta')) #glob of fasta files
			[shutil.copy(rf,sample.dirpath) for rf in reference_fasta_list] # copy to sample dir


			if no_deduplicate == False:
				# read deuplication. quality, adaptor trimming with fastp
				sample_logger.logger.info('Deduplicate reads, read trimming and quality scoring\n')
				call_Command(cmd=
					[
					'fastp', 
					'--in1', pjoin(sequenceDataDir,sample.read1_filename), #readfileR1
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
					'--thread', '3',
					'--dedup'
					],
					logger_=sample_logger)

			else :
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


			# for each reference available get coverage stats. produces the <reference>.fasta_coverage_stats.csv file
			sample_logger.logger.info('Calculate coverage\n')
			for reference in glob.glob('*.fasta'):
				calculate_ReferenceCoverage(sequenceDataDir=sequenceDataDir,reference=reference, read1_filename=sample.read1_filename, read2_filename=sample.read2_filename,samplename=sample.samplename, logger_=sample_logger)

			sample_logger.logger.info('Summarize coverage and pick reference\n')
			# summarize the results of strain selection. will be used in sample reports.
			summarize_ReferenceCoverage(samplename=sample.samplename)
			summarize_HACoverage(samplename=sample.samplename)
			combined_ReferenceCoverage(samplename=sample.samplename)
			# select the best reference strain and assign it to variables used by assemble.R
			reference_strain = select_BestReference(samplename=sample.samplename).replace('.fasta_coverage_stats.csv','.fasta')

			sample_logger.logger.info('Remove intermediate output files\n')
			# cleanup output files
			cleanup_CalculateReferenceCoverage(samplename=sample.samplename)

		except:
			sample_logger.logger.exception('FluPipeline Error: failure at strain assignment', exc_info=True)
			run_chunk = False

	
	# set new trimmed filenames
	trimmed_r1 = 'fastp_trimmed_'+sample.read1_filename
	trimmed_r2 = 'fastp_trimmed_'+sample.read2_filename


	## ==================================assemble reads==============================================================
	## ==============================================================================================================
	if run_chunk == True:
		try:
			if no_assembly == False:
				sample_logger.logger.info('Assembling reads')

				# run spades
				os.makedirs('assembly', exist_ok=True)
				call_Command(cmd=
					[
					'spades.py',
					'-1','{}'.format(pjoin(sample.dirpath,trimmed_r1)),
					'-2','{}'.format(pjoin(sample.dirpath,trimmed_r2)),
					'-o','{}/assembly'.format(sample.dirpath),
					'-t','4',
					'--phred-offset','33'
					],logger_=sample_logger
					)
			else:
				sample_logger.logger.info('Assembly option not selected')


		except:
			sample_logger.logger.exception('FluPipeline Error: failure at read assembly', exc_info=True)
			run_chunk = False			

	## ==================================find variants===============================================================
	## ==============================================================================================================
	if run_chunk == True:
		try:
			sample_logger.logger.info('Finding variants')

			# trimmed_r1 = 'fastp_trimmed_'+sample.read1_filename
			# trimmed_r2 = 'fastp_trimmed_'+sample.read2_filename

			# run assemble.R
			call_Command(cmd=
				[
				'{}'.format(Rscript),
				'{}/detect_variants.R'.format(softwareDir),
				'--workDir', '{}'.format(pjoin(sample.dirpath)),
				'--softwareDir', '{}'.format(softwareDir),
				'--R1', '{}'.format(trimmed_r1),
				'--R2', '{}'.format(trimmed_r2),
				'--refGenomeFasta', '{}'.format(reference_strain),
				'--minVariantPhredScore', '{}'.format(min_variant_phred_score),
				'--removeNTsFromAlignmentEnds', '{}'.format(remove_NTs_from_alignment_ends),
				'--BWAmappingScore', '{}'.format(min_read_mapping_score),
				'--minorVariantThreshold', '{}'.format(min_variant_frequency)
				],logger_=sample_logger
				)

		except:
			sample_logger.logger.exception('FluPipeline Error: failure at variant calling', exc_info=True)
			run_chunk = False

	if os.path.exists(detect_variants_message_file) == True:
		with open(detect_variants_message_file,'r') as infile:
			vardetect_message =  infile.readlines()[0].replace('\n','')
			sample_logger.logger.info('detect_variants.R message: '+vardetect_message)

		# these errors mean the rest of the chunks cannot be run
		if vardetect_message == 'No pileup was created':
			run_chunk = False
		if vardetect_message == 'No reads were aligned':
			run_chunk = False


	## ==================================Make consensus sequence=====================================================
	## ==============================================================================================================
	if run_chunk == True:

		try:
			sample_logger.logger.info('Making consensus sequences')
			# intialize consensus sequence object
			consensus = ConsensusSeq(samplename=sample.samplename, reference_strain=reference_strain)

			# check whether consensus sequence can be built from a variantTableMajor or not.
			if os.path.exists(sample.samplename+'_variantTableMajor.csv') == True:
				if len(pd.read_csv(sample.samplename+'_variantTableMajor.csv')) != 0:
					sample_logger.logger.info('variantTableMajor has variants')
					df_consensus = consensus.make_ConsensusTable(masking_threshold=consensus_masking_threshold)

				else:
					# if no variantTableMajor is available then the consensus sequence will be identical to the reference fasta except in the masked consensus sequence
					sample_logger.logger.info('variantTableMajor does not have variants')
					df_consensus = consensus.combine_PileupAndReferenceFasta()
					df_consensus['consensus']=df_consensus['REFERENCE']
					df_consensus['masked_consensus'] = np.where(df_consensus['pileup_cov']>consensus_masking_threshold,df_consensus['consensus'],'N')

				consensus.write_ConsensusFasta(df=df_consensus, base_from_col='consensus',filename_suffix='.consensus')
				consensus.write_ConsensusFasta(df=df_consensus, base_from_col='masked_consensus',filename_suffix='.consensus.masked')

				# store the consensus table that was used to make the consensus sequence
				df_consensus.to_csv(sample.samplename+'_consensusTable.csv',index=None)

		except:
			sample_logger.logger.exception('FluPipeline Error: failure at consensus sequence', exc_info=True)
			run_chunk = False


	## ==================================nextclade clade assignment==================================================
	## ==============================================================================================================
	if run_chunk == True:
		try:
			if masked_nextclade == True:
				nextclade_input_fasta = sample.samplename+'.consensus.masked.fasta'
			else:
				nextclade_input_fasta = sample.samplename+'.consensus.fasta'

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
					'--in-order', 
					'--input-fasta', '{}_HA.consensus.fasta'.format(sample.samplename),
					'--input-dataset', pjoin(os.getcwd(),nextclade_reference),
					'--output-tsv', 'nextclade_output/{}_clade_assignment.tsv'.format(sample.samplename), 
					'--output-tree', 'nextclade_output/{}.auspice.json'.format(sample.samplename),
					'--output-dir', 'nextclade_output/', 
					'--output-basename', '{}'.format(sample.samplename)
					],
					logger_=sample_logger)
			else:
				sample_logger.logger.info('No nextclade reference strain exists for this sample strain')

		except:
			sample_logger.logger.exception('FluPipeline Error: failure at nextclade clade assignment', exc_info=True)
			run_chunk = False

	
	## ==================================Intra-host variant detection================================================
	## ==============================================================================================================
	if run_chunk == True:
		try:
			sample_logger.logger.info('Finding intrahost variants')

			if masked_ivar == True:
				ivar_input_fasta = sample.samplename+'.consensus.masked.fasta'
			else:
				ivar_input_fasta = sample.samplename+'.consensus.fasta'

			# run assemble.R
			call_Command(cmd=
				[
				'{}'.format(Rscript),
				'{}/detect_variants.R'.format(softwareDir),
				'--workDir', '{}'.format(pjoin(intrahost_var_directory)),
				'--softwareDir', '{}'.format(softwareDir),
				'--R1', '{}'.format(pjoin(sample.dirpath, trimmed_r1)),
				'--R2', '{}'.format(pjoin(sample.dirpath, trimmed_r2)),
				'--refGenomeFasta', '{}'.format(pjoin(sample.dirpath,ivar_input_fasta)),
				'--minVariantPhredScore', '{}'.format(min_variant_phred_score),
				'--removeNTsFromAlignmentEnds', '{}'.format(remove_NTs_from_alignment_ends),
				'--BWAmappingScore', '{}'.format(min_read_mapping_score),
				'--minorVariantThreshold', '{}'.format(min_variant_frequency)
				],logger_=sample_logger
				)

		except:
			sample_logger.logger.exception('FluPipeline Error: failure at intrahost variant calling', exc_info=True)
			sample_logger.logger.info('FluPipeline Time: {}'.format(datetime.now()-start_run_timer))
			run_chunk = False

	if os.path.exists(pjoin(intrahost_var_directory,detect_variants_message_file)) == True:
		with open(pjoin(intrahost_var_directory,detect_variants_message_file),'r') as infile:
			vardetect_message =  infile.readlines()[0].replace('\n','')
			sample_logger.logger.info('detect_variants.R message: '+vardetect_message+' ivar')
			run_chunk = False


	## ==================================Sample report generation====================================================
	## ==============================================================================================================

	warnings = []
	with open(pjoin(baseDirectory,'sampleLogs',sample.samplename+'.txt'), 'r') as infile:
		for l in infile:
			if l.find('FluPipeline Error')>-1:
				print(sample.samplename,': ',l)
				warnings.append(l.replace('\n',''))
			if l.find('detect_variants.R message:')>-1:
				warnings.append(l.replace('\n',''))
	pd.DataFrame(warnings).to_csv('warnings.csv',index=None) #will write an empty file if there are no warnings


	try:
		sample_logger.logger.info('Generating sample report')
		# generate sample report
		call_Command(cmd=
			[
			'{}'.format(Rscript),
			'{}/report_runner.R'.format(softwareDir),
			'--softwareDir','{}'.format(softwareDir),
			'--report_type', '{}'.format('sample'), 
			'--sampleDir', '{}'.format(sample.dirpath), 
			'--samplename', '{}'.format(sample.samplename)
			],logger_=sample_logger
			)

	except:
		sample_logger.logger.exception('FluPipeline Error: failure at sample report generation', exc_info=True)

	# ==================================clean up and record run=====================================================
	# ==============================================================================================================
	
	if keep_all_intermediate_files == False:
		sample.cleanup_OutputFiles(directory=sample.dirpath)
		sample.cleanup_OutputFiles(directory=intrahost_var_directory)
		sample_logger.logger.info('Removed intermediate files')

	if run_chunk == True:
		sample_logger.logger.info('Finished processing sample {}'.format(sample.samplename))
	else:
		sample_logger.logger.info('Finished processing sample {} with warnings'.format(sample.samplename))

	sample_logger.logger.info('FluPipeline Time: {}'.format(datetime.now()-start_run_timer))
	#use print so the program talks more
	print('Finished processing sample {}'.format(sample.samplename))
	print('FluPipeline Time: {}'.format(datetime.now()-start_run_timer))

