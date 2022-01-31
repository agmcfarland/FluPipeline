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
global script_path
script_path = os.path.split(os.path.realpath(os.path.abspath(__file__)))[0]
sys.path.append(script_path)
from pylibs.processing_functions import *
from pylibs.processing_classes import SequencingSample, RunLogger, SampleTracker
pd.set_option('display.max_columns', None)

def run_FluPipeline(args):
	'''
	'''

	baseDirectory = args.base_directory
	sequenceDataDir = args.sequence_directory
	force_base_directory = args.force_base_directory
	referenceStrainsDir = pjoin(script_path,'references') #hardcoded default
	softwareDir = os.getcwd()
	threads = args.threads

	if args.cleanup == True:
		cleanup_files = args.cleanup
	else:
		cleanup_files = False

	if args.force == True:
		force_overwrite = args.force
	else:
		force_overwrite = False

	### INPUTS start ####'
	if force_base_directory == True:
		try:
			shutil.rmtree(baseDirectory)
		except:
			pass

	# input scripts
	Rscript = 'Rscript'

	#input other
	pipeline_used = 'snp' #options: snp, phylo
	process_method = 'bushman_artic_v2'
	#force_overwrite = True #deletes and remakes the sample folder in sampleOutputs
	BWA_path = subprocess.Popen(['which','bwa'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()[0].decode('ascii').replace('\n','')
	samtoolsbin_path = os.path.dirname(subprocess.Popen(['which','samtools'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()[0].decode('ascii').replace('\n',''))
	bcftoolsbin_path = os.path.dirname(subprocess.Popen(['which','bcftools'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()[0].decode('ascii').replace('\n',''))

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
			bcftoolsbin_path,
			cleanup_files
			])

	# log inputs
	run_logger.add_Message('inputs:') #contains all path and paramater inputs
	[run_logger.add_Message('{}: {}'.format(k, v)) for k,v in {
	'cleanup_files':cleanup_files, 'force_overwrite':force_overwrite,'force_base_directory':force_base_directory,
	'baseDirectory':baseDirectory,'Rscript':Rscript,'softwareDir':softwareDir,'pipeline_used':pipeline_used,
	'sequenceDataDir':sequenceDataDir,'referenceStrainsDir':referenceStrainsDir,'force_overwrite':force_overwrite,'BWA_path':BWA_path,'samtoolsbin_path':samtoolsbin_path
	}.items()]

	run_logger.add_Message('reference strains:')

	[run_logger.add_Message(os.path.basename(g)) for g in glob.glob(pjoin(referenceStrainsDir,'*.gb'))]

	run_logger.add_Message('samples:') # contains sample name, fastq R1 file, fastq R2 file
	for s in glob.glob(pjoin(sequenceDataDir,'*_R1_*.fastq.gz')):
		sample = SequencingSample()
		sample.get_DataFromReadPairs(read1_filename=s)
		run_logger.add_Message([sample.samplename,'  ',sample.read1_filename,'  ',sample.read2_filename])		

	run_logger.add_Message('total samples: {}'.format(len(sample_submission))) # total samples ran
	run_logger.add_Message('Processing samples...')

	## start data processing-----------------------------

	# convert genbank files to uniformally-formatted fasta. see convert_GBKTOFasta for details.
	gbk_reference_files = glob.glob(pjoin(referenceStrainsDir,'*.gb'))
	[convert_GBKToFasta(filename=f.replace('.gb','')) for f in gbk_reference_files]


	## run listed samples through specified worflow
	with Pool(processes=threads) as p:
		df_processed = pd.concat(p.starmap(flu_Pipeline, sample_submission))

	run_logger.add_Message('Finished processing samples...')
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



def main(args=None):
	'''
	Main function for parsing arguments
	'''

	if args is None:
		args = sys.argv[1:]

	# create parser
	# parser = argparse.ArgumentParser(prog='fluPipeline')
	parser = argparse.ArgumentParser()

	parser.add_argument('--base_directory',type=str ,default=None, help='directory that run samples will be saved in')
	parser.add_argument('--reference_directory',type=str ,default=None, help='directory containing reference strain files (.gb format)')
	parser.add_argument('--sequence_directory',type=str ,default=None, help='directory containing fastq sequence files (.gz format) ')
	parser.add_argument('--force', action='store_true', help='overwrite existing files')
	parser.add_argument('--force_base_directory', action='store_true', help='overwrite existing directory')
	parser.add_argument('--cleanup', action='store_true', help='remove intermediate files')
	parser.add_argument('--threads',type=int, default=4, help='number of processors to use for multiprocessing')
	parser.add_argument('--runtest', action='store_true', help='run an in silico test to make sure FluPipeline is working correctly')

	args = parser.parse_args()

	if args.runtest == True:
		print('yes')

	else:
		run_FluPipeline(args=args)

if __name__ == '__main__':
	sys.exit(main())


