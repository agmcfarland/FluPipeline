
import sys
import argparse
import logging
import pandas as pd
import shutil
import os
from os.path import join as pjoin
import subprocess
from Bio.Seq import Seq
from Bio import SeqIO
import glob
from pylibs.helper_functions import call_Command
from pylibs.processing_classes import RunLogger



if __name__ == '__main__':
	args = sys.argv[1:]
	parser = argparse.ArgumentParser(prog='best_reference')
	parser.add_argument('--baseDir', type=str, default=os.getcwd(), help='none', metavar='none')
	parser.add_argument('--logDir', type=str, default=os.getcwd(), help='none', metavar='none')
	parser.add_argument('--referenceStrainsDir', type=str, default=None, help='none', metavar='none')
	parser.add_argument('--R1', type=str, default=None, help='none', metavar='none')
	parser.add_argument('--R2', type=str, default=None, help='none', metavar='none')

	args = parser.parse_args()

	baseDir = args.baseDir
	logDir = args.logDir
	referenceStrainsDir = args.referenceStrainsDir
	R1 = args.R1
	R2 = args.R2

	# define sample name and specify explicit path
	os.chdir(baseDir)
	samplename = os.path.basename(os.getcwd())

	# start logging
	logger = RunLogger(directory=logDir,filename=samplename)
	logger.initialize_FileHandler()
	logger.logger.info('Running best_reference.py\n')
	logger.logger.info('Processing sample: {}'.format(samplename))


	logger.logger.info('Program arguments:\n')
	arguments_list = vars(args)
	for k,v in arguments_list.items():
		logger.logger.info('{}: {}\n'.format(k,v))

	# for each reference available get coverage stats. produces the <reference>.fasta_coverage_stats.csv file
	[shutil.copy(reffasta, os.getcwd()) for reffasta in glob.glob(referenceStrainsDir+'/*.fasta')]

	# check for important data
	if os.path.exists(R1) == False:
		raise ValueError(f'R1 does not exist:]\n{R1}')
	if os.path.exists(R2) == False:
		raise ValueError(f'R2 does not exist:\n{R2}')


	# calculate coverage for each reference file
	for reference in glob.glob(baseDir+'/*.fasta'):

		reference = os.path.basename(reference)

		## index reference
		call_Command(cmd=
			[
			'bwa','index',reference
			],
			logger_=logger)
		# break
		## align random read subset to reference with bwa
		call_Command(cmd=
			'bwa mem -M {} {} {} > {}_bwa_alignment_{}.sam'.format(reference, R1, R2, reference, samplename)
			,
			logger_=logger,
			shell_=True)

		## sam to bam and sort and index alignment using samtools
		call_Command(cmd=
			'samtools view -S -b {}_bwa_alignment_{}.sam > {}_bwa_alignment_{}.bam'.format(reference,samplename,reference,samplename),
			logger_=logger,
			shell_=True)
		call_Command(cmd=
			'samtools sort -o {}_bwa_alignment_{}.bam {}_bwa_alignment_{}.bam'.format(reference,samplename,reference,samplename),
			logger_=logger,
			shell_=True)
		call_Command(cmd=
			'samtools index {}_bwa_alignment_{}.bam'.format(reference,samplename),
			logger_=logger,
			shell_=True)

		## pileup data with samtools
		call_Command(cmd=
			'samtools mpileup -A -a -Q 0 -o {}_pileup_{}.txt -d 100000 -f {} {}_bwa_alignment_{}.bam'.format(reference,samplename, reference, reference,samplename),
			logger_=logger,
			shell_=True)


		## calculate coverage per segment and average depth per segment
		df = pd.read_csv('{}_pileup_{}.txt'.format(reference,samplename), names=['segment','position','base','coverage','quality','extra'], engine='python', sep='\t', quoting=3)

		if len(df) > 0:

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

		else:
			df_out = pd.DataFrame({
				'segment': [1,2,3,4,5,6,7,8],
				'segment_coverage': [0,0,0,0,0,0,0,0],
				'average_read_depth': [0,0,0,0,0,0,0,0],
				})

		df_out['reference'] = reference # add the reference used for tracking purposes

		## store mapping results
		df_out.to_csv('{}_coverage_stats.csv'.format(reference),index=False)


	# Aggregate summarized coverage statistics
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

		df_temp = pd.DataFrame({
			'reference':[reference_strain],
			'segment_count':[segment_count],
			'segment_coverage_average':[segment_coverage_average],
			'segment_average_read_depth':[segment_average_read_depth],
			})

		if len(df_cov) == 0:
			df_cov = df_temp
		else:
			df_cov = pd.concat([df_cov,df_temp])


	df_cov = df_cov.sort_values(['segment_count','segment_coverage_average','segment_average_read_depth'], ascending = [False,False,False])
	df_cov['sample'] = samplename
	df_cov.to_csv('reference_coverage_{}.csv'.format(samplename), index=False)


	# combine coverage stats into one long dataframe
	df_cov = pd.DataFrame()
	for reference_strain in glob.glob('*_coverage_stats.csv'):
		df = pd.read_csv(reference_strain)
		if len(df_cov) == 0:
			df_cov = df
		else:
			df_cov = pd.concat([df_cov, df])
	df_cov.to_csv('combined_coverage_stats_{}.csv'.format(samplename), index=False)

	# summarize HA coverage stats
	df_ha = df_cov[df_cov['segment']==4]
	df_ha = df_ha.sort_values(['segment_coverage','average_read_depth'], ascending = [False,False])
	df_ha['sample'] = samplename
	df_ha.to_csv('reference_ha_coverage_{}.csv'.format(samplename), index=False)

	# cleanup files
	df = pd.read_csv('reference_ha_coverage_{}.csv'.format(samplename))

	# the reference strains are sorted with the best match at the top. select the the first item in the reference column
	reference_to_keep = df['reference'][0].replace('.fasta_coverage_stats.csv','')


	for f in os.listdir():
		if f.find('fastp') > -1:
			continue
		if f.find('reference_ha_coverage') > -1:
			continue
		if f.find('reference_coverage') > -1:
			continue
		if f.find('combined_coverage_stats') > -1:
			continue
		if os.path.isdir(f) == True:
			continue
		else:
			os.remove(f)






