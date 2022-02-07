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


class SequencingSample:
	'''
	Treats each sequencing sample as an object with a sample name, two read files, and an directory containing all data created
	In
	Additional methods for listing the read files and checking their existence
	'''
	def __init__(self):
		'''
		Initiates the object
		'''
		self.read1_filename = ''
		self.read2_filename = ''
		self.samplename = ''
		self.dirpath = ''

	def get_DataFromReadPairs(self, read1_filename):
		'''
		Uses the R1 read to assign r2 filename and extract the sample name
		'''
		self.read1_filename = os.path.basename(read1_filename)
		self.read2_filename = self.read1_filename.replace('_R1_','_R2_')
		self.samplename = re.findall(r'.*_R1',self.read1_filename)[0].replace('_R1','')

	def get_DataFromSampleName(self, samplename, readfile_suffix='_R1_001.fastq.gz'):
		'''
		Uses the samplename to assign r1 and r2 filenames. Readfile_suffix is optional
		'''
		self.samplename = samplename
		self.read1_filename = self.samplename+readfile_suffix
		self.read2_filename = self.samplename+readfile_suffix

	def get_SampleDirPath(self, directory):
		'''
		Assigns the sample directory path using a supplied base directory.
		Default is the current directory/sampleOutputs.
		'''
		self.dirpath = pjoin(directory,'sampleOutputs', self.samplename)

	def check_ReadExistence(self, directory):
		'''
		Checks that the reads in the object exist in the specified directory
		'''
		if os.path.exists(pjoin(directory,self.read1_filename)) == False:
			raise ValueError(str(self.read1_filename)+' not found')
		if os.path.exists(pjoin(directory,self.read2_filename)) == False:
			raise ValueError(str(self.read2_filename)+' not found')

	def cleanup_OutputFiles(self, clean_level='intermediate'):
		'''
		Removes files after running pipeline that do not have a partial string match to intermediate_keep list
		intermediate: all intermediate files 
		none: no files removed
		'''
		intermediate_keep = [
		'.bam','filt.qual.sorted.bam.bai','filt.qual.sorted.bam','vcf.gz','.filt.vcf.gz','.filt.vcf.gz.csi',
		'pileup','.Rdata','.pdf','fastp_stats','fastp.html','reference_coverage','_coverage_stats.csv','variantTable',
		'reference_ha_coverage','consensus.fasta'
		]

		to_remove = os.listdir(self.dirpath)

		for f in os.listdir(self.dirpath):
			for ik in intermediate_keep:
				if f.find(ik) > -1:
					try:
						to_remove.remove(f)
					except:
						continue

		if clean_level == 'intermediate':
			[os.remove(f) for f in to_remove]



class RunLogger:
	'''
	Handles logging: Making a unique logfile by name. Adding messages to it. Removing logfiles.
	'''
	def __init__(self, directory, filename):
		'''
		directory is where the logfile will be stored
		filename is the name of the file. no suffix is supplied
		'''
		self.directory = directory
		self.filename = filename
		self.filepath = pjoin(self.directory,self.filename)
		self.formatter = logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(message)s',
			datefmt='%Y-%m-%d %H:%M:%S')

	def initialize_FileHandler(self):
		'''
		Initialize a logfile with time stamps and in append mode.
		'''
		handler = logging.FileHandler('{}.txt'.format(self.filepath), mode='a')
		handler.setFormatter(self.formatter)
		self.logger = logging.getLogger(self.filename)
		self.logger.setLevel(logging.DEBUG)
		self.logger.addHandler(handler)

	def add_StreamHandler(self):
		'''
		Adds the stream handler to the logger object
		'''
		handler = logging.StreamHandler()
		handler.setFormatter(self.formatter)
		self.logger.addHandler(handler)


if __name__=='__main__':
	pass









