import sys
import os 
from os.path import join as pjoin
import multiprocessing as mp
from multiprocessing import Pool
import pandas as pd
import numpy as np
import re
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

	def cleanup_OutputFiles(self, directory):
		'''
		Removes files after running pipeline that do not have a partial string match to intermediate_keep list
		intermediate: all intermediate files 
		none: no files removed
		'''
		curdir = os.getcwd()
		os.chdir(directory)

		keep = [
			'.consensus.fasta','consensusTable.csv','consensus.masked.fasta','variantTableRaw.csv',
			'variantTableMajor.csv','variantTable.csv','variantQualityHisto.txt','allVariants.vcf',
			'covstats.txt','basecov.txt','filt.qual.sorted.bam.bai','filt.qual.sorted.bam',
			'reference_ha_coverage','reference_coverage','combined_coverage_stats',
			'coverage_stats','fastp_stats',
			'nextclade_output','_ha','detect_variants_messages.txt','warnings.csv','.pdf'
		]
		# change so it doesn't remove folders
		to_remove = os.listdir(directory)

		for f in os.listdir(directory):
			for i in keep:
				if f.find(i) > -1:
					try:
						to_remove.remove(f)
					except:
						continue

		# [print(f) for f in to_remove if os.path.isdir(f)==False]
		[os.remove(f) for f in to_remove if os.path.isdir(f)==False]

		os.chdir(curdir)

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



class ConsensusSeq:
	'''
	Creates a consensus sequence using pileup and variant tables
	'''

	def __init__(self, samplename, reference_strain):
		self.samplename = samplename
		self.reference_strain = reference_strain
		if reference_strain.endswith('.fasta')==False:
			raise ValueError('reference strain must have a .fasta extension')

	def make_ReferenceTable(self):
		'''
		Takes as input a fasta file.
		Makes a dataframe of the reference bases for each segment and their nucleotide position (starts at one for each segment)
		'''
		df_reference = pd.DataFrame()
		for record in SeqIO.parse(self.reference_strain,'fasta'):
			df_reference = df_reference.append([[record.id,i] for i in str(record.seq)])
		df_reference.columns = ['CHROM','REFERENCE']
		df_reference['POS'] = df_reference.groupby('CHROM').cumcount()+1 #+1 so that counting starts at 1
		return(df_reference)

	def format_Pileup(self):
		'''
		Takes as input a basecov pileup file
		Makes adataframe with the per segment base positions
		'''
		df_pup = pd.read_csv(pjoin(self.samplename+'_basecov.txt'),sep='\t')
		df_pup.columns = ['CHROM','POS','pileup_cov']
		df_pup['POS'] = df_pup['POS']+1 #+1 so that counting starts at 1
		return(df_pup)

	def combine_PileupAndReferenceFasta(self):
		'''
		Merges the outputs from make_ReferenceTable and format_Pileup
		'''
		df_reference = self.make_ReferenceTable()
		df_pup = self.format_Pileup()

		# Check lengths are equivalent. If not, then something went wrong during the coverage calculation process.
		if len(df_pup)!=len(df_reference):
			raise ValueError('Pileup coverage file and fasta file must have same segment lengths.')

		# Merge the reference with the consensus to get coverage and reference bases
		df_combined = df_reference.merge(df_pup,on=['CHROM','POS'])

		return(df_combined)

	def make_ConsensusTable(self, masking_threshold=0):
		'''
		Takes as input a consensus table from combine_PileupAndReferenceFasta
		Returns a table of CHROM, REFERENCE, POS, pileup_voc, consensus, masked consensus
		'''
		if len(pd.read_csv(pjoin(self.samplename+'_variantTableMajor.csv'))) == 0:
		# if len(pd.read_csv(pjoin(samplename+'_variantTableMajor.csv'))) == 0:	
			raise ValueError('no variants in variantTableMajor.csv')

		# Filter for deletions in the variantTableMajor file and make a long-form datatable for each deletion range (i.e. 22-25 becomes 23,24,25)
		df_mv = pd.read_csv(pjoin(self.samplename+'_variantTableMajor.csv'))[['CHROM','POS','REF','ALT','TYP','STA','STO']]
		# df_mv = pd.read_csv(pjoin(samplename+'_variantTableMajor.csv'))[['CHROM','POS','REF','ALT','TYP','STA','STO']]
		df_del = pd.DataFrame()
		for index, row in df_mv.iterrows():
			if row['TYP'] == 'DEL':
				for i in range(row['STA']+1,row['STO']+1): # add 1 to each start/stop value. The deletion positions are listed after the start
					df_del = df_del.append([[row['CHROM'],row['TYP'],i]])

		# Merge consensus dataframe with deletion dataframe using an outer join. Then filter for bases that are on the left. This removes all bases that are deleted according the df_del
		df_consensus = self.combine_PileupAndReferenceFasta()

		if len(df_del) != 0:
			df_del.columns = ['CHROM','TYP','POS']
			df_consensus = df_consensus.merge(df_del,on=['CHROM','POS'],how='outer',indicator=True)
			df_consensus = df_consensus[df_consensus['_merge']=='left_only']
			df_consensus = df_consensus.drop(['_merge','TYP'],axis=1) # remove _merge and TYP 

		# Make a dataframe of insertion or substitution variants
		df_vars = df_mv[(df_mv['TYP']=='INS')|(df_mv['TYP']=='SUB')][['CHROM','POS','REF','ALT']]

		# Merge with consensus using outer join
		df_consensus = df_consensus.merge(df_vars, on=['CHROM','POS'],how='outer',indicator=True)

		# The consensus column will use the ALT (variant) nucelotides if df_consensus and df_vars had an entry for that segment&nucleotide. Else, it'll use the original reference nucleotide.
		df_consensus['consensus'] = np.where(df_consensus['_merge']=='both',df_consensus['ALT'],df_consensus['REFERENCE'])

		# The masked consensus uses the consensus base except in instances where the pile_up coverage is less than an inputted value
		df_consensus['masked_consensus'] = np.where(df_consensus['pileup_cov']>masking_threshold,df_consensus['consensus'],'N')

		df_consensus = df_consensus.drop(['_merge','REF','ALT'],axis=1) #remove unneeded columns for final output

		return(df_consensus)

	def write_ConsensusFasta(self, df, base_from_col, filename_suffix):
		'''
		Takes as input either a df_consensus output from ConsensusSeq.make_ConsensusTable or
		df_combined output from ConsensusSeq.combine_PileupAndReference to write a new table
		'''
		with open(self.samplename+filename_suffix+'.fasta','w') as outfile:
			for chrom in df['CHROM'].unique().tolist():
				df_temp = df[df['CHROM']==chrom]
				chrom_seq = ''.join(df_temp[base_from_col].tolist())
				outfile.write('>{}\n{}\n'.format(chrom, chrom_seq))



if __name__=='__main__':
	pass









