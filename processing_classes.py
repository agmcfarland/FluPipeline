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
		'reference_ha_coverage'
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
	Handles logging.
	'''
	def __init__(self, directory, filename):
		'''
		directory is where the logfile will be stored
		filename is the name of the file. no suffix is supplied
		'''
		self.directory = directory
		self.filename = filename
		self.filepath = pjoin(self.directory,self.filename)


	def initialize_Log(self):
		'''
		Initialize a logfile with time stamps and in append mode.
		'''
		formatter = logging.Formatter(fmt='%(asctime)s %(levelname)-8s %(message)s',
			datefmt='%Y-%m-%d %H:%M:%S')

		handler = logging.FileHandler('{}.txt'.format(self.filepath), mode='a')
		handler.setFormatter(formatter)
		self.logger = logging.getLogger(self.filename)
		self.logger.setLevel(logging.DEBUG)
		self.logger.addHandler(handler)


	def add_Message(self, message, level='info'):
		'''
		Add a message to the logfile
		'''
		if level == 'info':
			self.logger.info(message)
		if level == 'debug':
			self.logger.debug(message)
		if level == 'warning':
			self.logger.warning(message)

	def delete_LogFile(self):
		'''
		deletes the logfile
		'''
		if os.path.exists('{}.txt'.format(self.filepath)):
			os.remove('{}.txt'.format(self.filepath))


class SampleTracker:
	'''
	Loads the sample tracker sql database into the object. Check existence of tables, pull tables, update tables in server, save tables in .csv format.
	'''
	def __init__(self, directory):
		'''
		'''
		self.conn = sql.connect(pjoin(directory,'sampleTracker.db'))

	def check_TableExists(self, tablename):
		'''
		Check that a given processing table exists. if not then create it with the specified columns
		'''
		cursor = self.conn.cursor()
		try:
			cursor.execute('SELECT sample from {}'.format(tablename))
			print('table {} exists'.format(tablename))
		except:
			cursor.execute("""CREATE TABLE {} (sample, processed, date_run, runtime, info)""".format(tablename))
			print('created table {}'.format(tablename))
		cursor.close()

	def load_Table(self, tablename):
		'''
		Load a table with a given name into a pandas dataframe
		'''
		return(pd.read_sql_query('SELECT * FROM {}'.format(tablename), self.conn))

	def update_Table(self, tablename, df, if_exists_):
		'''
		Update a sql table with a given name by supplying it a pandas dataframe.
		'''
		if tablename not in ['snp','phylo']:
			raise ValueError('choose either snp or phylo')
		df.to_sql(name=tablename,con=self.conn, if_exists=if_exists_,index_label='sample',index=False)

	def write_Table(self, tablename, name, output_dir=os.getcwd()):
		'''
		Writes a table in the database with a given name to a file in a specificied directory
		'''
		df = pd.read_sql_query('SELECT * FROM {}'.format(tablename), self.conn)
		df.to_csv(pjoin(output_dir,name+'.csv'),index=None)











