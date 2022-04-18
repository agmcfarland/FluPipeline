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
import psutil

global script_path
script_path = os.path.split(os.path.realpath(os.path.abspath(__file__)))[0]
sys.path.append(script_path)
from pylibs.processing_functions import *
from pylibs.processing_classes import SequencingSample, RunLogger, ConsensusSeq
from pylibs.pipeline import flu_Pipeline
from pylibs.insilico_test import create_TestData
from pylibs.helper_functions import automatic_ThreadUsage, call_Command
pd.set_option('display.max_columns', None)

def run_FluPipeline(args):
	'''
	Runs through all steps of FluPipeline
	'''
	start_run_timer = datetime.now()
	## prepare working directory
	if args.sequence_directory == pjoin(script_path,'references'): # prevent the user from using the default reference directory in combination with the --use_fasta flag
		args.use_fasta = False

	if args.force_base_directory == True: 	# remove base directory if it exiss
		try:
			shutil.rmtree(args.base_directory)
		except:
			pass

	## make inputs for assemble.R
	# assigning softwareDir as a copy of script_path for logical consistency with assemble.R inputs.
	Rscript = 'Rscript'

	## set up main directories
	os.makedirs(args.base_directory, exist_ok=True)
	os.chdir(args.base_directory)
	os.makedirs(pjoin(args.base_directory,'sampleLogs'), exist_ok=True)
	os.makedirs(pjoin(args.base_directory,'sampleOutputs'), exist_ok=True)
	os.makedirs(pjoin(args.base_directory,'sampleReports'), exist_ok=True)

	## start up run logger
	run_logger = RunLogger(directory=os.getcwd(),filename='runLog')
	run_logger.initialize_FileHandler()
	run_logger.add_StreamHandler()
	run_logger.logger.info('Starting FluPipeline...\n')

	## log computing resources available and requested
	run_logger.logger.info('Total memory available: {}Gb\n'.format(round(psutil.virtual_memory()[0]/(1024.0**3)),2))
	run_logger.logger.info('Threads available: {}\n'.format(os.cpu_count()))
	if args.max_mem_per_thread != None: # check whether automatic thread number was requested and log calculations
		run_logger.logger.info('Automatic thread number selection\n')
		run_logger.logger.info('Threads use: {}\n'.format(args.threads))
		run_logger.logger.info('Estimated memory per thread: {}Gb\n'.format(args.max_mem_per_thread))
		run_logger.logger.info('Estimated max memory usage: {}Gb\n'.format(args.threads*args.max_mem_per_thread))

	## log argument inputs, reference strains, and samples	
	run_logger.logger.info('Program arguments:\n')
	arguments_list = vars(args)
	for k,v in arguments_list.items():
		run_logger.logger.info('{}: {}\n'.format(k,v))

	run_logger.logger.info('Reference strains:\n')
	if args.use_fasta == True:
		[run_logger.logger.info('{}\n'.format(os.path.basename(g))) for g in glob.glob(pjoin(args.reference_directory,'*.fasta'))]
	else:
		[run_logger.logger.info('{}\n'.format(os.path.basename(g))) for g in glob.glob(pjoin(args.reference_directory,'*.gb'))]

	run_logger.logger.info('Samples:\n') # contains sample name, fastq R1 file, fastq R2 file
	for s in glob.glob(pjoin(args.sequence_directory,'*_R1_*.fastq.gz')):
		sample = SequencingSample()
		sample.get_DataFromReadPairs(read1_filename=s)
		run_logger.logger.info('{}\n'.format(sample.samplename))


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
			args.keep_all_intermediate_files,
			args.strain_sample_depth,
			args.downsample,
			args.consensus_masking_threshold,
			args.min_variant_phred_score,
			args.remove_NTs_from_alignment_ends,
			args.min_read_mapping_score,
			args.masked_nextclade,
			args.masked_ivar,
			args.base_quality,
			args.no_deduplicate			
			])

	run_logger.logger.info('Total samples to process: {}\n'.format(len(sample_submission))) # total samples ran

	if len(sample_submission) == 0:
		run_logger.logger.exception('No samples available to process. Check directory for correctly formatted fastq file names.\n')

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


	## write software versions
	run_logger.logger.info('Recording software/package versions...\n')
	call_Command(cmd=
		'conda list --export > softwareVersions.txt',
		logger_=run_logger,
		shell_=True)
		

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

	# copy sample reports pdfs to sampleReports directory
	[shutil.copy(report, pjoin(args.base_directory,'sampleReports')) for report in glob.glob(pjoin(args.base_directory,'sampleOutputs','*','*.pdf'))]


	## end run
	run_logger.logger.info('Finished running FluPipeline\n')
	run_logger.logger.info('Total time: {}\n'.format(datetime.now()-start_run_timer))
	run_logger.logger.info('Run ouputs stored in {}\n'.format(args.base_directory))	



def main(args=None):
	'''
	Main function for parsing arguments
	'''
	if args is None:
		args = sys.argv[1:]

	# create parser
	parser = argparse.ArgumentParser(prog = 'FluPipeline')

	# add arguments to parser
	parser.add_argument('--base_directory',type=str ,default=None, help='directory that run samples will be saved in. [None]', metavar='')
	parser.add_argument('--reference_directory',type=str ,default=pjoin(script_path,'references'), help='directory containing reference strain files (.gb or .fasta (see --use_fasta flag)) [script_path/references]', metavar='')
	parser.add_argument('--sequence_directory',type=str ,default=None, help='directory containing fastq sequence files (.gz format) [None]', metavar='')
	parser.add_argument('--force', action='store_true', default=False, help='overwrite existing sample files. [False]')
	parser.add_argument('--force_base_directory', action='store_true', default=False, help='overwrite existing directory. [False')
	parser.add_argument('--max_mem_per_thread', type=int, default=None, help='automatically determines the number of threads to use based on memory per thread supplied (in Gb) [None]', metavar='')
	parser.add_argument('--keep_all_intermediate_files', action='store_true', default=False, help='remove intermediate files. [False]')
	parser.add_argument('--threads',type=int, default=4, help='number of samples to process in paralell. one sample is one read pair [4]', metavar='')
	parser.add_argument('--strain_sample_depth', type=int, default=2000, help='number of random reads to use to determine strain assignment. [2000]', metavar='')
	parser.add_argument('--use_fasta', action='store_true', default=False, help='fasta format: fasta file(s) contain all eight segments sequences. All segments must have a single name (only letters, numbers, and underscores. At the end of the name there should be an underscore followed by the segment number. Example: an_example_name_1. [False]')
	parser.add_argument('--consensus_masking_threshold', type=int, default=0, help='replace any nucleotides in the consensus sequence with N if their depth falls below this number. [0]', metavar='')
	parser.add_argument('--downsample', type=int, default=-1, help='downsample all read files to these many reads. [-1 (no downsampling)]', metavar='')
	parser.add_argument('--min_variant_phred_score', type=int, default=5, help='keep all variants above or equal to this phred-scaled value. [5]', metavar='')
	parser.add_argument('--remove_NTs_from_alignment_ends', type=int, default=3, help='remove this many bases from the left and right of each read prior to mapping. [3]', metavar='')
	parser.add_argument('--min_read_mapping_score', type=int, default=30, help='keep reads that mapped above or eequal to this phred-scaled value. [3]', metavar='')
	parser.add_argument('--masked_nextclade', action='store_true', default=False, help='use the masked consensus sequence fasta file for nextclade clade assignment.  [False]')
	parser.add_argument('--masked_ivar', action='store_true', default=False, help='use the masked consensus sequence fasta file for intrahost variation detection.  [False]')
	parser.add_argument('--base_quality', type=int, default=30, help='keep reads that have at least an average of this pphred-scaled value. [30]', metavar='')
	parser.add_argument('--no_deduplicate',  action='store_true', default=False, help='do not conduct read deduplication.  [False]')
	parser.add_argument('--runtest', action='store_true', default=False, help='run an in silico test to make sure FluPipeline is working correctly. [False]')


	# create args opbject with arguments
	args = parser.parse_args()

	## calculate threads to use if --max_mem_per_thread was supplied
	if args.max_mem_per_thread != None:
		calculated_threads_to_use = automatic_ThreadUsage(max_mem_per_thread = args.max_mem_per_thread)
		args.threads = calculated_threads_to_use #overwrite args.threads with the new value

	## create test data if --runtest argument was supplied
	if args.runtest == True:
		testDir =  pjoin(script_path,'run_test') #change script path to os.getcwd() in the future.
		args.base_directory = pjoin(testDir,'output')
		args.sequence_directory = pjoin(testDir,'data')
		args.reference_directory = pjoin(script_path,'references')
		args.force == True
		args.force_base_directory = False
		args.keep_all_intermediate_files = True
		create_TestData(testDir=testDir, referenceStrainsDir=args.reference_directory)

	##---Run FluPipline ---##	
	run_FluPipeline(args=args)

if __name__ == '__main__':
	sys.exit(main())

