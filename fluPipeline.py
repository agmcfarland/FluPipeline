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
from pylibs.processing_classes import SequencingSample, RunLogger
from pylibs.pipeline import flu_Pipeline
from pylibs.insilico_test import *
pd.set_option('display.max_columns', None)

def run_FluPipeline(args):
	'''
	Runs through all steps of FluPipeline
	'''
	## prepare working directory
	# prevent the user from using the default reference directory in combination with the use_fasta
	if args.sequence_directory == pjoin(script_path,'references'):
		args.use_fasta = False

	# remove base directory if it exiss
	if args.force_base_directory == True:
		try:
			shutil.rmtree(args.base_directory)
		except:
			pass


	## make inputs for assemble.R
	# assigning softwareDir as a copy of script_path for logical consistency with assemble.R inputs.
	Rscript = 'Rscript'
	process_method = 'bushman_artic_v2'
	BWA_path = subprocess.Popen(['which','bwa'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()[0].decode('ascii').replace('\n','')
	samtoolsbin_path = os.path.dirname(subprocess.Popen(['which','samtools'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()[0].decode('ascii').replace('\n',''))
	bcftoolsbin_path = os.path.dirname(subprocess.Popen(['which','bcftools'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()[0].decode('ascii').replace('\n',''))


	## set up main directories
	os.makedirs(args.base_directory, exist_ok=True)
	os.chdir(args.base_directory)
	os.makedirs(pjoin(args.base_directory,'sampleLogs'), exist_ok=True)
	os.makedirs(pjoin(args.base_directory,'sampleOutputs'), exist_ok=True)
	os.makedirs(pjoin(args.base_directory,'sampleResults'), exist_ok=True)

	## start up run logger
	run_logger = RunLogger(directory=os.getcwd(),filename='runLog')
	run_logger.initialize_FileHandler()
	run_logger.add_StreamHandler()
	run_logger.logger.info('Starting FluPipeline...\n')

	## log argument inputs, reference strains, and samples	
	run_logger.logger.info('inputs:\n') #contains all path and paramater inputs
	run_logger.logger.info(args)

	run_logger.logger.info('reference strains:\n')
	if args.use_fasta == True:
		[run_logger.logger.info(os.path.basename(g)) for g in glob.glob(pjoin(args.reference_directory,'*.fasta'))]
	else:
		[run_logger.logger.info(os.path.basename(g)) for g in glob.glob(pjoin(args.reference_directory,'*.gb'))]

	run_logger.logger.info('samples:\n') # contains sample name, fastq R1 file, fastq R2 file
	for s in glob.glob(pjoin(args.sequence_directory,'*_R1_*.fastq.gz')):
		sample = SequencingSample()
		sample.get_DataFromReadPairs(read1_filename=s)
		run_logger.logger.info(sample.samplename)		


	## start data processing------------------------------------------------------------

	## collect all samples to process in a list for multiprocessing submission
	sample_submission = []
	for f in glob.glob(pjoin(args.sequence_directory,'*_R1_*.fastq.gz')):
		# create sample object with necessary data and check that it exists
		sample = SequencingSample()
		sample.get_DataFromReadPairs(read1_filename=f)
		sample.get_SampleDirPath(directory=args.base_directory)
		sample.check_ReadExistence(directory=args.sequence_directory)

		# append sample run parameters to master run list
		sample_submission.append([
			args.base_directory,
			Rscript,
			script_path,
			sample,
			args.sequence_directory,
			args.reference_directory,
			args.force,
			BWA_path,
			samtoolsbin_path,
			bcftoolsbin_path,
			args.cleanup,
			args.strain_sample_depth
			])

	run_logger.logger.info('Total samples to process: {}\n'.format(len(sample_submission))) # total samples ran

	if len(sample_submission) == 0:
		run_logger.logger.exception('No samples available to process. Check directory for correctly formatted fastq file names.\n')
		sys.exit('Stopped FluPipeline')

	if args.use_fasta == True:
		run_logger.logger.info('Using fasta files for reference\n')
	else:
		run_logger.logger.info('Using gbk files for reference\n')
		# convert genbank files to uniformally-formatted fasta. see convert_GBKTOFasta for details.
		gbk_reference_files = glob.glob(pjoin(args.reference_directory,'*.gb'))
		[convert_GBKToFasta(filename=f.replace('.gb','')) for f in gbk_reference_files]


	run_logger.logger.info('Processing samples...\n')
	## run listed samples through flu_Pipeline function
	with Pool(processes=args.threads) as p:
		p.starmap(flu_Pipeline, sample_submission)

	run_logger.logger.info('Finished processing samples\n')

	## end data processing--------------------------------------------------------------


	## store run metadata
	run_logger.logger.info('Storing run statistics in runStats.csv\n')
	# gather run metadata from logfile in a list and then store as a pandas df
	run_metadata = []
	for s in glob.glob(pjoin(args.sequence_directory,'*_R1_*.fastq.gz')):
		sample = SequencingSample()
		sample.get_DataFromReadPairs(read1_filename=s)
		with open(pjoin(os.getcwd(),'sampleLogs',sample.samplename+'.txt'),'r') as infile:
			error_status = ''
			time_ran = ''
			for l in infile:
				if l.find('FluPipeline Error:') > -1:
					error_status = 'error'
				if l.find('FluPipeline Time:') > -1:
					time_ran = l[l.find('me: ')+4:].replace('\n','') #use 'me: ' as an anchor + 4 positions to extract just the time
			run_metadata.append([datetime.now(),sample.samplename,error_status,time_ran,pjoin(args.sequence_directory,sample.read1_filename),pjoin(args.sequence_directory,sample.read2_filename)])
	df_meta = pd.DataFrame(run_metadata)
	df_meta.columns = ['date_time','sample','error_status','time_ran','read1_filepath','read2_filepath']

	# save the run metadata to file. if there already is a file, the metadata get appeneded to it.
	if os.path.exists('runStats.csv') == True:
		df_meta1 = pd.read_csv('runStats.csv')
		df_meta1 = df_meta1.append(df_meta)
		df_meta1.to_csv('runStats.csv',index=None)
	else:
		df_meta.to_csv('runStats.csv',index=None)


	## make run report
	run_logger.logger.info('Making run summary...\n')
	call_Command(cmd=
		[
		'{}'.format(Rscript),
		 '{}/report_runner.R'.format(script_path),
		 '--softwareDir', '{}'.format(script_path),
		 '--report_type', '{}'.format('run'),
		 '--baseDir', '{}'.format(args.base_directory)
		 ], logger_=run_logger)

	run_logger.logger.info('Finished making run summary\n')

	# copy sample reports pdfs to sampleResults directory
	[shutil.copy(report, pjoin(args.base_directory,'sampleResults')) for report in glob.glob(pjoin(args.base_directory,'sampleOutputs','*','*.pdf'))]


	## end run
	run_logger.logger.info('Run ouputs stored in {}\n'.format(args.base_directory))	
	run_logger.logger.info('Finished running FluPipeline\n')



def main(args=None):
	'''
	Main function for parsing arguments
	'''

	if args is None:
		args = sys.argv[1:]

	# create parser
	parser = argparse.ArgumentParser()

	# add arguments to parser
	parser.add_argument('--base_directory',type=str ,default=None, help='directory that run samples will be saved in')
	parser.add_argument('--reference_directory',type=str ,default=pjoin(script_path,'references'), help='directory containing reference strain files (.gb or .fasta (see --use_fasta flag))')
	parser.add_argument('--sequence_directory',type=str ,default=None, help='directory containing fastq sequence files (.gz format) ')
	parser.add_argument('--force', action='store_true', default=False, help='overwrite existing files in assemble.R script')
	parser.add_argument('--force_base_directory', action='store_true', default=False, help='overwrite existing directory')
	parser.add_argument('--cleanup', action='store_true', default=False, help='remove intermediate files')
	parser.add_argument('--threads',type=int, default=4, help='number of processors to use for multiprocessing')
	parser.add_argument('--runtest', action='store_true', default=False, help='run an in silico test to make sure FluPipeline is working correctly')
	parser.add_argument('--strain_sample_depth', type=int, default=2000, help='number of random reads to use to determine strain assignment. default=2000')
	parser.add_argument('--use_fasta', action='store_true', default=False, help='Fast format: fasta file(s) contain all eight segments sequences. All segments must have a single name (only letters, numbers, and underscores. At the end of the name there should be an underscore followed by the segment number. Example: an_example_name_1. default=False')

	# create args opbject with arguments
	args = parser.parse_args()

	##---Run test #---
	if args.runtest == True:
		testDir =  pjoin(script_path,'run_test') #change script path to os.getcwd() in the future.
		args.base_directory = pjoin(testDir,'output')
		args.sequence_directory = pjoin(testDir,'data')
		args.reference_directory = pjoin(script_path,'references')
		args.force == True
		args.force_base_directory = False
		args.cleanup = False
		# args.threads = 6
		create_TestData(testDir=testDir, referenceStrainsDir=args.reference_directory)
		run_FluPipeline(args=args)

	##---Run FluPipleine #---	
	else:
		run_FluPipeline(args=args)

if __name__ == '__main__':
	sys.exit(main())


