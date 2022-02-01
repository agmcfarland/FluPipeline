
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
from .processing_classes import SampleTracker



def create_SyntheticReads(filename):
	'''
	crete synthetic reads from fasta file using bbmap
	'''
	# # # no errors
	# subprocess.run(['randomreads.sh',
	# 	'reads=100000',
	# 	# 'coverage=1000',
	# 	'prefix={}'.format(filename),
	# 	'seed=5','paired=t','addpairnum=f','illuminanames=t','adderrors=f',
	# 	'ref={}.fasta'.format(filename),'out1={}_perfect_R1_001.fastq'.format(filename),'out2={}_perfect_R2_001.fastq'.format(filename),
	# 	'minlength=150','maxlength=150',
	# 	'maxsnps=0','snprate=0',
	# 	'maxinss=0',
	# 	'maxdels=0',
	# 	'maxsubs=0',
	# 	'q=40',
	# 	'superflat=t',
	# 	'overwrite=t'])
	
	# # SNPs
	subprocess.run(['randomreads.sh',
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
		'overwrite=t'])

	# # # insertions deletions and SNPs
	# subprocess.run(['randomreads.sh',
	# 	'reads=100000',
	# 	# 'coverage=1000',
	# 	'prefix={}'.format(filename),
	# 	'seed=5','paired=t','addpairnum=f','illuminanames=t','adderrors=f',
	# 	'ref={}.fasta'.format(filename),'out1={}_snpindel_R1_001.fastq'.format(filename),'out2={}_snpindel_R2_001.fastq'.format(filename),
	# 	'minlength=150','maxlength=150',
	# 	'maxsnps=2','snprate=1',
	# 	'maxinss=2','insrate=1','maxinslen=3',
	# 	'maxdels=2','delrate=1','maxdellen=3',
	# 	'q=40',
	# 	'superflat=t',
	# 	'overwrite=t'])



def generate_TestData(filename):
	'''
	Supply the filename with no suffix (i.e. remove .fasta)
	'''
	convert_GBKToFasta(filename=filename)
	create_SyntheticReads(filename=filename)



def create_TestData(testDir, referenceStrainsDir):
	'''
	Main function for this script
	'''
	# testDir = '/home/agmcfarland/flu_project/FluPipeline_test/'
	baseDirectory = pjoin(testDir,'output')
	sequenceDataDir = pjoin(testDir,'data')
	# referenceStrainsDir = '/home/agmcfarland/flu_project/FluPipeline/references'


	# remove testing directory if it already exists
	try:
		shutil.rmtree(testDir)
	except:
		pass

	# make testing directories
	os.makedirs(testDir, exist_ok=True)
	os.makedirs(baseDirectory, exist_ok=True)
	os.makedirs(sequenceDataDir, exist_ok=True)

	# enter the test directory that will contain the insilico fastq files
	os.chdir(sequenceDataDir)

	# copy reference genbank files (.gb) from referenceStrainsDir to the sequence data dir, sequenceDataDir
	[shutil.copy(gbk, os.getcwd()) for gbk in glob.glob(pjoin(referenceStrainsDir,'*.gb'))]

	# extract just the filename, no suffix for all gbs to convert to synthetic reads
	gbk_files = glob.glob('*.gb')
	gbk_files = [i.replace('.gb','') for i in gbk_files]

	# create fasta and synthetic read set
	for strain in gbk_files:
		generate_TestData(filename=strain)

	# check perfect, snp, and indel are all made
	print(len(gbk_files) == len(glob.glob('*perfect*'))/2)
	print(len(gbk_files) == len(glob.glob('*_snp_*'))/2)
	print(len(gbk_files) == len(glob.glob('*snpindel*'))/2)

	# gzip fastq files
	trimmed_fastq_list = glob.glob('*.fastq')
	[subprocess.run(['gzip', trimmed_fastq]) for trimmed_fastq in trimmed_fastq_list] 


## script end
if __name__=='__main__':
	pass







