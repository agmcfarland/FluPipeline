#!/usr/bin/env python

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
pd.set_option('display.max_columns', None)


global script_path
script_path = os.path.split(os.path.realpath(os.path.abspath(__file__)))[0]
sys.path.append(script_path)


from processing_functions import *
from processing_classes import SequencingSample, RunLogger, SampleTracker


def run_FluPipeline(args):
	'''
	'''

	### INPUTS start ####
	try:
		shutil.rmtree('/home/agmcfarland/flu_project/test/test4')
	except:
		pass
	# input data
	# sequenceDataDir = '/home/agmcfarland/flu_project/test/test_data' #'/home/ashley/Fastq' # sequencing data comes from this directory
	sequenceDataDir = '/home/agmcfarland/flu_project/test/test_data' #'/home/ashley/Fastq' # sequencing data comes from this directory
	baseDirectory = '/home/agmcfarland/flu_project/test/test4' # all data/directory tree will be saved in this directory
	referenceStrainsDir = '/home/agmcfarland/flu_project/test/test_data' #directory reference strains to find the best match for a given readset

	# input scripts
	Rscript = 'Rscript'#'/home/opt/R-3.4.0/bin/Rscript'
	softwareDir = '/home/agmcfarland/flu_project/flu_pipeline'

	#input other
	pipeline_used = 'snp' #options: snp, phylo
	process_method = 'bushman_artic_v2'
	force_overwrite = True #deletes and remakes the sample folder in sampleOutputs
	BWA_path =  '/home/agmcfarland/miniconda3/envs/FluPipeline_env/bin/bwa'#'/home/everett/ext/bwa' #to keep bwa version consistent
	samtoolsbin_path = '/home/agmcfarland/miniconda3/envs/FluPipeline_env/bin'#'/home/everett/ext/samtools/bin' #to keep samtools verison consistent
	bcftoolsbin_path = '/home/agmcfarland/miniconda3/envs/FluPipeline_env/bin'

	### INPUTS end ####

	# set up main directories
	os.makedirs(baseDirectory, exist_ok=True)
	os.chdir(baseDirectory)
	os.makedirs(pjoin(baseDirectory,'sampleLogs'), exist_ok=True)
	os.makedirs(pjoin(baseDirectory,'sampleOutputs'), exist_ok=True)
	os.makedirs(pjoin(baseDirectory,'sampleResults'), exist_ok=True)

	# start up run logger
	run_logger = RunLogger(directory=os.getcwd(),filename='runLog')
	run_logger.delete_LogFile() # deletes a previously existing logfile with the same name
	run_logger.initialize_Log()
	run_logger.add_Message('Start run')

	## load/create sample tracker sql database and tables for tracking
	sample_tracker = SampleTracker(directory=baseDirectory)
	sample_tracker.check_TableExists(tablename=pipeline_used)
	df_st = sample_tracker.load_Table(tablename=pipeline_used) #load sql table into pandas df 

	print(df_st)
	## collect all samples to process in a list for multiprocessing submission
	sample_submission = []
	for f in glob.glob(pjoin(sequenceDataDir,'*_R1_*.fastq.gz')):
		# create sample object with necessary data and check that it exists
		sample = SequencingSample()
		sample.get_DataFromReadPairs(read1_filename=f)
		sample.get_SampleDirPath(directory=baseDirectory)
		sample.check_ReadExistence(directory=sequenceDataDir)

		# if sample.samplename.find('perfect') == -1:
		# 	continue

		# append sample run parameters to master run list
		sample_submission.append([
			baseDirectory,
			Rscript,
			softwareDir,
			sample,
			pipeline_used,
			sequenceDataDir,
			referenceStrainsDir,
			force_overwrite,
			BWA_path,
			samtoolsbin_path,
			bcftoolsbin_path
			])

	# log inputs
	run_logger.add_Message('inputs:') #contains all path and paramater inputs
	[run_logger.add_Message('{}: {}'.format(k, v)) for k,v in {'baseDirectory':baseDirectory,'Rscript':Rscript,'softwareDir':softwareDir,'pipeline_used':pipeline_used,'sequenceDataDir':sequenceDataDir,'referenceStrainsDir':referenceStrainsDir,'force_overwrite':force_overwrite,'BWA_path':BWA_path,'samtoolsbin_path':samtoolsbin_path}.items()]
	run_logger.add_Message('reference strains:')
	[run_logger.add_Message(os.path.basename(g)) for g in glob.glob(pjoin(referenceStrainsDir,'*.gb'))]
	run_logger.add_Message('samples:') # contains sample name, fastq R1 file, fastq R2 file
	for s in glob.glob(pjoin(sequenceDataDir,'*_R1_*.fastq.gz')):
		sample = SequencingSample()
		sample.get_DataFromReadPairs(read1_filename=s)
		run_logger.add_Message([sample.samplename,'  ',sample.read1_filename,'  ',sample.read2_filename])
	run_logger.add_Message('total samples: {}'.format(len(sample_submission))) # total samples ran


	## start data processing-----------------------------

	# convert genbank files to uniformally-formatted fasta. see convert_GBKTOFasta for details.
	gbk_reference_files = glob.glob(pjoin(referenceStrainsDir,'*.gb'))
	[convert_GBKToFasta(filename=f.replace('.gb','')) for f in gbk_reference_files]


	## run listed samples through specified worflow
	with Pool(processes=10) as p:
		df_processed = pd.concat(p.starmap(flu_Pipeline, sample_submission))


	print('time to make report!!!!!')
	run_logger.add_Message('Making run report...')
	## make run report
	os.system('{} {}/report_runner.R --softwareDir {} --report_type {} --baseDir {}'.
		format(
			Rscript,
			softwareDir,
			softwareDir, #--softwareDir
			'run', #--report_type
			baseDirectory #--baseDir
			))
	run_logger.add_Message('Finished making run report')

	# gather sample reports
	[shutil.copy(report, pjoin(baseDirectory,'sampleResults')) for report in glob.glob(pjoin(baseDirectory,'sampleOutputs','*','*.pdf'))]

	## end data processing-------------------------------

	## update the sample tracker database
	sample_tracker.update_Table(tablename=pipeline_used, df=df_processed, if_exists_='append')
	sample_tracker.write_Table(tablename=pipeline_used, name='runStats', output_dir=baseDirectory)

	## close connection to sql database
	sample_tracker.conn.close()

	run_logger.add_Message('Finished run')



def main(args = None):
	'''
	Main function for parsing arguments
	'''

	if args is None:
		args = sys.argv[1:]

	print(args)
	# create parser
	# parser = argparse.ArgumentParser(prog='fluPipeline')
	parser = argparse.ArgumentParser()

	parser.add_argument('--base_directory',type=str ,default=None, help='directory that run samples will be saved in')
	parser.add_argument('--reference_directory',type=str ,default=None, help='directory containing reference strain files (.gb format)')
	parser.add_argument('--sequence_directory',type=str ,default=None, help='directory containing fastq sequence files (.gz format) ')
	parser.add_argument('--force',type=int, default=False, help='overwrite existing files')
	parser.add_argument('--threads',type=int, default=None, help='number of processors to use for multiprocessing')

	args = parser.parse_args()

	run_FluPipeline(args=args)

if __name__ == '__main__':
	sys.exit(main())


