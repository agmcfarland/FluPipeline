import sys
import os 
from os.path import join as pjoin
import multiprocessing as mp
from multiprocessing import Pool
import pandas as pd
import re
import sqlite3 as sql
import glob
import shutil
import subprocess
from datetime import datetime
from Bio.Seq import Seq
from Bio import SeqIO, SeqFeature
import argparse
import logging
import json
from .processing_classes import SequencingSample, RunLogger
from .processing_functions import *

def flu_Pipeline(
	baseDirectory,
	Rscript,
	softwareDir,
	sample,
	sequenceDataDir,
	referenceStrainsDir,
	force_overwrite,
	BWA_path,
	samtoolsbin_path,
	bcftoolsbin_path,
	cleanup_files,
	strain_sample_depth):
	'''
	Runs through all steps of the flu pipeline. 
	'''

# ## troubleshooting inputs start ##
# # use this exact input for troubleshooting the pipeline
# baseDirectory='/home/agmcfarland/flu_project/FluPipeline/run_test/output'
# Rscript='Rscript'
# softwareDir = '/home/agmcfarland/flu_project/FluPipeline'
# sample = SequencingSample()
# sample.get_DataFromReadPairs(read1_filename='IBV_Yamagata_ref_snp_R1_001.fastq.gz')
# sample.get_SampleDirPath(directory=baseDirectory)
# pipeline_used='snp'
# sequenceDataDir = '/home/agmcfarland/flu_project/FluPipeline/run_test/data'
# referenceStrainsDir = '/home/agmcfarland/flu_project/FluPipeline/references' #directory reference strains to find the best match for a given readset
# pipeline_used = 'snp' #options: snp, phylo
# process_method = 'bushman_artic_v2'
# force_overwrite = True #deletes and remakes the sample folder in sampleOutputs
# BWA_path =  '/home/agmcfarland/miniconda3/envs/FluPipeline_env/bin/bwa'#'/home/everett/ext/bwa' #to keep bwa version consistent
# samtoolsbin_path = '/home/agmcfarland/miniconda3/envs/FluPipeline_env/bin'#'/home/everett/ext/samtools/bin' #to keep samtools verison consistent
# bcftoolsbin_path = '/home/agmcfarland/miniconda3/envs/FluPipeline_env/bin'
# cleanup_files = True
# strain_sample_depth = 3000


# use this exact input for troubleshooting the pipeline
# baseDirectory='/home/agmcfarland/flu_project/shared_data/test_6_samples'
# Rscript='Rscript'
# softwareDir = '/home/agmcfarland/flu_project/FluPipeline'
# sample = SequencingSample()
# sample.get_DataFromReadPairs(read1_filename='ashley_5_R1_001.fastq.gz')
# sample.get_SampleDirPath(directory=baseDirectory)
# pipeline_used='snp'
# sequenceDataDir = '/home/agmcfarland/flu_project/shared_data/test_data_6_samples'
# referenceStrainsDir = '/home/agmcfarland/flu_project/FluPipeline/references' #directory reference strains to find the best match for a given readset
# pipeline_used = 'snp' #options: snp, phylo
# process_method = 'bushman_artic_v2'
# force_overwrite = True #deletes and remakes the sample folder in sampleOutputs
# BWA_path =  '/home/agmcfarland/miniconda3/envs/FluPipeline_env/bin/bwa'#'/home/everett/ext/bwa' #to keep bwa version consistent
# samtoolsbin_path = '/home/agmcfarland/miniconda3/envs/FluPipeline_env/bin'#'/home/everett/ext/samtools/bin' #to keep samtools verison consistent
# bcftoolsbin_path = '/home/agmcfarland/miniconda3/envs/FluPipeline_env/bin'
# cleanup_files = True
# strain_sample_depth = 3000
# ## troubleshooting inputs end ##


	## ==================================Prepare sample run directory================================================
	## ==============================================================================================================
	# set and create directory to store all intermediate data in sample-specific folders. Switch to that directory afterwards.

	if force_overwrite==True:
		try: 
			shutil.rmtree(sample.dirpath)
		except:
			pass

	os.makedirs(sample.dirpath, exist_ok=True)
	os.chdir(sample.dirpath)

	start_run_timer = datetime.now()
	sample_logger = RunLogger(directory=pjoin(baseDirectory,'sampleLogs'),filename=sample.samplename)
	sample_logger.initialize_FileHandler()
	sample_logger.logger.info('Starting run')
	sample_logger.logger.info('sample: {}'.format(sample.samplename))

	## ==================================Find and assign reference strain============================================
	## ==============================================================================================================
	try:
		sample_logger.logger.info('Assigning reference strain')
		# copy reference fastas to sample directory 
		reference_fasta_list = glob.glob(pjoin(referenceStrainsDir,'*.fasta')) #glob of fasta files
		[shutil.copy(rf,sample.dirpath) for rf in reference_fasta_list] # copy to sample dir


		# quality, adaptor trimming with fastp
		sample_logger.logger.info('Step 1: read trimming and quality scoring\n')
		call_Command(cmd=
			[
			'fastp', 
			'--in1', pjoin(sequenceDataDir,sample.read1_filename), #readfileR1
			'--in2', pjoin(sequenceDataDir,sample.read2_filename), #readfileR2
			'--out1', 'fastp_trimmed_{}'.format(sample.read1_filename), 
			'--out2', 'fastp_trimmed_{}'.format(sample.read2_filename),
			'-j', 'fastp_stats_{}'.format(sample.samplename), #json output
			'-q', '30',#, #quality score 30
			'overwrite=True'#,
			#'--disable_adapter_trimming' # adaptor trimming when adaptors are already trimmed can lead to unwanted errors
			],
			logger_=sample_logger)

		# convert fastp json into table
		with open('fastp_stats_{}'.format(sample.samplename)) as infile:
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
		sample_logger.logger.info('Step 2: subsample reads\n')
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
		sample_logger.logger.info('Step 3: calculate coverage\n')
		for reference in glob.glob('*.fasta'):
			calculate_ReferenceCoverage(sequenceDataDir=sequenceDataDir,reference=reference, BWA_path=BWA_path, samtoolsbin_path=samtoolsbin_path, read1_filename=sample.read1_filename, read2_filename=sample.read2_filename,samplename=sample.samplename, logger_=sample_logger)

		sample_logger.logger.info('Step 4: summarize coverage and pick reference\n')
		# summarize the results of strain selection. will be used in sample reports.
		summarize_ReferenceCoverage(samplename=sample.samplename)
		summarize_HACoverage(samplename=sample.samplename)
		combined_ReferenceCoverage(samplename=sample.samplename)
		# select the best reference strain and assign it to variables used by assemble.R
		reference_strain = select_BestReference(samplename=sample.samplename).replace('.fasta_coverage_stats.csv','.fasta')
		refGenomeFasta = reference_strain
		refGenomeBWA = reference_strain

		sample_logger.logger.info('Step 5: remove intermediate output files\n')
		# cleanup output files
		cleanup_CalculateReferenceCoverage(samplename=sample.samplename)

	except:
		sample_logger.logger.exception('FluPipeline Error: failure at strain assignment', exc_info=True)
		sample_logger.logger.info('FluPipeline Time: {}'.format(datetime.now()-start_run_timer))

	## ==================================assemble.R find variants====================================================
	## ==============================================================================================================
	try:
		sample_logger.logger.info('Finding variants')

		# run assemble.R
		call_Command(cmd=
			[
			'{}'.format(Rscript),
			'{}/assemble.R'.format(softwareDir),
			'--outputFile', '{}'.format(pjoin(sample.dirpath,sample.samplename+'.Rdata')),
			'--workDir', '{}'.format(pjoin(sample.dirpath)),
			'--process_method', '{}'.format('bushman_artic_v2'),
			'--softwareDir', '{}'.format(softwareDir),
			'--R1', '{}'.format(pjoin(sequenceDataDir,sample.read1_filename)),
			'--R2', '{}'.format(pjoin(sequenceDataDir,sample.read2_filename)),
			'--refGenomeFasta', '{}'.format(refGenomeFasta),
			'--refGenomeBWA', '{}'.format(refGenomeBWA),
			'--bwaPath', '{}'.format(BWA_path),
			'--samtoolsBin', '{}'.format(samtoolsbin_path),
			'--bcftoolsBin', '{}'.format(bcftoolsbin_path)
			],logger_=sample_logger
			)

	except:
		sample_logger.logger.exception('FluPipeline Error: failure at variant calling', exc_info=True)
		sample_logger.logger.info('FluPipeline Time: {}'.format(datetime.now()-start_run_timer))

	
	## ==================================nextclade clade assignment==================================================
	## ==============================================================================================================
	try:
		nextclade_reference = reference_NextCladeLookUp(reference=reference_strain.replace('.fasta',''))

		if nextclade_reference != 'none available':
			sample_logger.logger.info('Finding clade assignment using nextclade')

			# make new fasta file containing only the HA sequence (segment 4) that is used by nextclade for clade assignment
			with open('{}_HA.consensus.fasta'.format(sample.samplename),'w') as infile:
				for record in SeqIO.parse('{}.consensus.fasta'.format(sample.samplename),'fasta'):
					if record.id.endswith('_4') == True:
						infile.write('>{}-{}\n'.format(sample.samplename,record.id))
						infile.write('{}\n'.format(str(record.seq)))

			# copy the reference directory files from the nextclade_references directory and store it in the sample directory			
			shutil.copytree(src=pjoin(softwareDir,'nextclade_references',nextclade_reference), dst=pjoin(os.getcwd(),nextclade_reference))

			# run nextclade. all ouputs are stored in a folder called 'nextclade_output'
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
		sample_logger.logger.info('FluPipeline Time: {}'.format(datetime.now()-start_run_timer))

	## ==================================Sample report generation====================================================
	## ==============================================================================================================

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
		sample_logger.logger.info('FluPipeline Time: {}'.format(datetime.now()-start_run_timer))

	# ==================================clean up and record run=====================================================
	# ==============================================================================================================
	if cleanup_files == True:
		sample.cleanup_OutputFiles(clean_level='intermediate')
		sample_logger.logger.info('Removed intermediate files')

	sample_logger.logger.info('Finished processing sample')
	sample_logger.logger.info('FluPipeline Time: {}'.format(datetime.now()-start_run_timer))

