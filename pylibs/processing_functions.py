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
from .processing_classes import SequencingSample, RunLogger


def call_Command(cmd, logger_, shell_=False):
	'''
	Runs a shell command using subprocess.run. Recods output to a logger object that has already been created.
	Default is to use subprocess.run. If shell is necessary then subprocess.Popen is used.
	'''
	try:
		if shell_ == False:
			capture = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True).stdout
			logger_.logger.info(capture)
		else:
			capture = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			logger_.logger.info(capture.communicate()[0].decode('utf-8'))
	except:
		logger_.logger.exception('\nFluPipeline Error:', exc_info=True)


def convert_GBKToFasta(filename):
	'''
	Creates a fasta file from genbank that uses the accession ids and segment numbers to name each segment. 
	'''
	with open('{}.fasta'.format(filename),'w') as infile:

		for record in SeqIO.parse('{}.gb'.format(filename),'gb'):

			record_id = record.id.replace(' ','_').replace('.','_')

			for f in record.features:
				if f.type == 'source':
					segment_id = f.qualifiers['segment'][0] 

			infile.write('>{}\n'.format(record_id+'_'+segment_id))
			infile.write('{}\n'.format(str(record.seq)))


def calculate_ReferenceCoverage(sequenceDataDir, reference, BWA_path, samtoolsbin_path, read1_filename, read2_filename, samplename, logger_):
	'''
	Runs through all steps of a simple read mapping pipeline and calculates the coverage and depth for a given read file reference combo.
	Ouputs a csv containing the coverage stats for all segments
	'''

# ## troubleshooting inputs start ##
# read1_filename = sample.read1_filename
# read2_filename = sample.read2_filename
# samplename = sample.samplename
# reference = 'H1N1pdm_ref.fasta'
# ## troubleshooting inputs end ##

	## index reference
	call_Command(cmd=
		[
		BWA_path,'index',reference
		],
		logger_=logger_)

	## align random read subset to reference with bwa
	call_Command(cmd=
		'{} mem -M {} subset_fastp_trimmed_{} subset_fastp_trimmed_{} > {}_bwa_alignment_{}.sam'.format(BWA_path, reference, read1_filename, read2_filename, reference, samplename)
		,
		logger_=logger_,
		shell_=True)

	## sam to bam and sort and index alignment using samtools
	call_Command(cmd=
		'{}/samtools view -S -b {}_bwa_alignment_{}.sam > {}_bwa_alignment_{}.bam'.format(samtoolsbin_path,reference,samplename,reference,samplename),
		logger_=logger_,
		shell_=True)
	call_Command(cmd=
		'{}/samtools sort -o {}_bwa_alignment_{}.bam {}_bwa_alignment_{}.bam'.format(samtoolsbin_path,reference,samplename,reference,samplename),
		logger_=logger_,
		shell_=True)
	call_Command(cmd=
		'{}/samtools index {}_bwa_alignment_{}.bam'.format(samtoolsbin_path,reference,samplename),
		logger_=logger_,
		shell_=True)

	## pileup data with samtools
	call_Command(cmd=
		'{}/samtools mpileup -A -a -Q 0 -o {}_pileup_{}.txt -d 100000 -f {} {}_bwa_alignment_{}.bam'.format(samtoolsbin_path, reference,samplename, reference, reference,samplename),
		logger_=logger_,
		shell_=True)


	## calculate coverage per segment and average depth per segment
	df = pd.read_csv('{}_pileup_{}.txt'.format(reference,samplename), names=['segment','position','base','coverage','quality','extra'], engine='python', sep='\t', quoting=3)

	df['segment'] = df['segment'].str.strip().str[-1] #rename segments to just their numbers 
	 
	df['binary_coverage'] = [0 if x == 0 else 1 for x in df['coverage'].tolist()] #any depth higher than 0 is 1

	df_seglengths = df[['segment','coverage']].groupby('segment').count().reset_index().rename(columns={'coverage':'segment_length'}) # genomic length of each segment

	df_segcoverage = df[['segment','binary_coverage']].groupby('segment').sum() # number of basepairs with a depth greater than 0

	df_averagecoverage = df[['segment','coverage']].groupby('segment').sum() # total reads at all positions

	# merge summary data into master dataframe df_out for calculations
	df_out = df_seglengths.merge(df_segcoverage, left_on='segment',right_on='segment') 

	df_out = df_out.merge(df_averagecoverage, left_on='segment',right_on='segment')

	df_out['segment_coverage'] = df_out['binary_coverage']/df_out['segment_length'] # percentage of bases with a coverage of 1 or greater per segment

	df_out['average_read_depth'] = df_out['coverage']/df_out['segment_length'] # average read depth per base per segment

	df_out = df_out[['segment','segment_coverage','average_read_depth']]

	df_out['reference'] = reference # add the reference used for tracking purposes

	## store mapping results
	df_out.to_csv('{}_coverage_stats.csv'.format(reference),index=False)


def summarize_ReferenceCoverage(samplename):
	'''
	Write summary files of coverage, depth, segment count, for all references compared. 
	Gives averages over all segments detected.
	'''
	df_cov = pd.DataFrame()
	for reference_strain in glob.glob('*_coverage_stats.csv'):
		df = pd.read_csv(reference_strain)

		if len(df) == 0:
			segment_count = 0
			segment_coverage_average = 0
			segment_average_read_depth = 0
		else:
			segment_count = len(df)
			segment_coverage_average = df['segment_coverage'].mean()
			segment_average_read_depth = df['average_read_depth'].mean()

		df_cov = df_cov.append(pd.DataFrame({
			'reference':[reference_strain],
			'segment_count':[segment_count],
			'segment_coverage_average':[segment_coverage_average],
			'segment_average_read_depth':[segment_average_read_depth],
			}), ignore_index = True)


	df_cov = df_cov.sort_values(['segment_count','segment_coverage_average','segment_average_read_depth'], ascending = [False,False,False])
	df_cov['sample'] = samplename
	df_cov.to_csv('reference_coverage_{}.csv'.format(samplename), index=False)


def combined_ReferenceCoverage(samplename):
	'''
	Writes a dataframe of combined *_coverage_stats.csv dataframes
	'''
	df_cov = pd.DataFrame()
	for reference_strain in glob.glob('*_coverage_stats.csv'):
		df = pd.read_csv(reference_strain)
		df_cov = df_cov.append(df)
	df_cov.to_csv('combined_coverage_stats_{}.csv'.format(samplename), index=False)



def summarize_HACoverage(samplename):
	'''
	Write summary files of coverage, depth, segment count, for all references compared. 
	'''
	df_ha = pd.DataFrame()
	for reference_strain in glob.glob('*_coverage_stats.csv'):
		df = pd.read_csv(reference_strain)

		if len(df[df['segment']==4]) > 0: #in case some segments have coverage but not at segment 4
			ha_coverage = df[df['segment']==4]['segment_coverage'].item()
			ha_read_depth = df[df['segment']==4]['average_read_depth'].item()
		else:
			ha_coverage = 0
			ha_read_depth = 0

		df_ha = df_ha.append(pd.DataFrame(
			{'reference':[reference_strain],
			'ha_coverage':[ha_coverage],
			'ha_read_depth':[ha_read_depth]
			}))

	df_ha = df_ha.sort_values(['ha_coverage','ha_read_depth'], ascending = [False,False])
	df_ha['sample'] = samplename
	df_ha.to_csv('reference_ha_coverage_{}.csv'.format(samplename), index=False)



def select_BestReference(samplename):
	'''
	Selects the best reference from the output of summarize_ReferenceCoverage
	(most segments, highest coverage, highest read depth). 
	'''
	df = pd.read_csv(pjoin('reference_ha_coverage_{}.csv'.format(samplename)))
	df = df.sort_values(['ha_coverage','ha_read_depth'], ascending = [False,False])
	return(df.head(1)['reference'].item())



def cleanup_CalculateReferenceCoverage(samplename):
	'''
	Removes all files except the reference summary csv files, the reference strain .fasta file, the reference strain BWA index,
	and the fastp statistics.

	Only reference strain index files, fastp summary stats, and reference coverage sumamry ,
	and each individual <reference>_coverage_stats.csv file are kept.
	'''
	df = pd.read_csv('reference_ha_coverage_{}.csv'.format(samplename))

	# the reference strains are sorted with the best match at the top. select the the first item in the reference column
	reference_to_keep = df['reference'][0].replace('.fasta_coverage_stats.csv','')


	# remove files that are not the reference
	for f in df['reference'][1:]:
		f = f.replace('.fasta_coverage_stats.csv','')
		for i in glob.glob('{}*'.format(f)):
			os.remove(i)
	# remove subsetted fastq and fastp trimmed fastq files
	for i in glob.glob('fastp_trimmed*'):
		os.remove(i)
	for i in glob.glob('subset_fastp*'):
		os.remove(i)
	# remove all sam and bam and bai files
	for i in glob.glob('*.bam'):
		os.remove(i)
	for i in glob.glob('*.bam.bai'):
		os.remove(i)
	for i in glob.glob('*.sam'):
		os.remove(i)


def reference_NextCladeLookUp(reference):
	'''
	Returns the the nextclade folder name containing all necessary comparison data to run nextclade
	'''
	reference_dict = {
	'IBV_Victoria_ref':'flu_vic_ha',
	'IBV_Yamagata_ref':'flu_yam_ha',
	'H1N1pdm_ref':'flu_h1n1pdm_ha',
	'H3N2_ref':'flu_h3n2_ha',
	'H6N2_Massachusetts':'none available',
	'H1N1_A_Brisbane_59_2007':'none available'
	}
	return(reference_dict[reference])



def flag_PotentialReassortment():
	'''
	'''
	pass



if __name__=='__main__':
	pass


