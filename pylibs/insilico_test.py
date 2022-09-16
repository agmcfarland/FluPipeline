
import subprocess
import os 
from os.path import join as pjoin
import pandas as pd
import re
from Bio.Seq import Seq
from Bio import SeqIO, SeqFeature
import glob
import shutil
import numpy as np
import logging
import multiprocessing as mp
from multiprocessing import Pool
from .processing_functions import convert_GBKToFasta
from .processing_classes import RunLogger
from .helper_functions import call_Command



def create_SyntheticReads(filename, logger_, shell_=False):
	'''
	crete synthetic reads from fasta file using bbmap
	'''
	call_Command(
		[
		'randomreads.sh',
		'reads=100000',
		# 'coverage=1000',
		'prefix={}'.format(filename),
		'seed=5','paired=t','addpairnum=f','illuminanames=t','adderrors=f',
		'ref={}.fasta'.format(filename),'out1={}_snp_R1_001.fastq'.format(filename),'out2={}_snp_R2_001.fastq'.format(filename),
		'minlength=150','maxlength=150',
		'maxsnps=2','snprate=1',
		'maxinss=0',
		'maxdels=0',
		'maxsubs=0',
		'q=40',
		'superflat=t',
		'overwrite=t'
		],logger_=logger_)


def create_TestData(testDir, referenceStrainsDir):
	'''
	Main function for this script
	'''
	baseDirectory = pjoin(testDir,'output')
	sequenceDataDir = pjoin(testDir,'data')

	# # remove testing directory if it already exists
	# try:
	# 	shutil.rmtree(testDir)
	# except:
	# 	pass

	# # make testing directories
	# os.makedirs(testDir, exist_ok=True)
	# os.makedirs(baseDirectory, exist_ok=True)
	# os.makedirs(sequenceDataDir, exist_ok=True)

	# # enter the test directory that will contain the insilico fastq files
	# os.chdir(sequenceDataDir)

	# # intialize logfile handler
	# run_logger = RunLogger(directory=baseDirectory,filename='insilicoTestLog')
	# run_logger.initialize_FileHandler()
	# run_logger.add_StreamHandler()

	# # copy reference genbank files (.gb) from referenceStrainsDir to the sequence data dir, sequenceDataDir
	# [shutil.copy(gbk, os.getcwd()) for gbk in glob.glob(pjoin(referenceStrainsDir,'*.gb'))]

	# # extract just the filename, no suffix for all gbs to convert to synthetic reads
	# gbk_files = glob.glob('*.gb')
	# gbk_files = [i.replace('.gb','') for i in gbk_files]

	# # create fasta and synthetic read set
	# run_logger.logger.info('Making test data...\n')
	# for strain in gbk_files:
	# 	convert_GBKToFasta(filename=strain)
	# 	create_SyntheticReads(filename=strain, logger_=run_logger)

	# # gzip fastq files
	# run_logger.logger.info('Compressing in silico fastq files...\n')
	# trimmed_fastq_list = glob.glob('*.fastq')
	# [subprocess.run(['gzip', trimmed_fastq]) for trimmed_fastq in trimmed_fastq_list] 
	# run_logger.logger.info('Finished making test data\n')


def compare_TestResults(): 
	'''
	Compares the results of the test routine with those produced by the development version of FluPipeline
	'''
	pass


## script end
if __name__=='__main__':
	pass







