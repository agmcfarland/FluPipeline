
import pandas as pd
import numpy as np
import psutil
from os.path import join as pjoin
import sys
import argparse
import logging
import shutil
import os
from os.path import join as pjoin
import subprocess
from Bio.Seq import Seq
from Bio import SeqIO
import glob
import functools
import json
from pylibs.processing_classes import SequencingSample
from pylibs.processing_functions import select_BestReference

pd.set_option('display.max_columns', None)


def safe_ReadCSV(filepath, sep_=','):
	'''
	'''
	if os.path.exists(filepath) == True:
		return pd.read_csv(filepath, sep = sep_)
	else:
		return pd.DataFrame()

def format_VariantTable(df):
	'''
	'''
	df = df[['CHROM','TYPE']].groupby(['CHROM','TYPE']).size().reset_index(name='COUNT')
	df = df.pivot(index='CHROM', columns='TYPE', values='COUNT')
	df['order'] = [i[-1] for i in list(df.index)]
	df = df.sort_values('order',ascending=True).drop(columns='order')
	return df

if __name__ == '__main__':

	args = sys.argv[1:]
	parser = argparse.ArgumentParser(prog='summarize')
	parser.add_argument('--baseDir', type=str, default=os.getcwd(), help='none', metavar='none')
	parser.add_argument('--sequenceDir', type=str, default=None, help='none', metavar='none')
	parser.add_argument('--use_strain', type=str, default=None, help='none', metavar='none')

	args = parser.parse_args()

	if args.use_strain == 'None':
		args.use_strain = None

	baseDir = args.baseDir
	sequenceDir = args.sequenceDir
	use_strain = args.use_strain

	# baseDir = '/data/Pipeline_output/033122_Output/033122_Output_LM'
	# sequenceDir = '/data/Fastq/033122_Fastq/033122-Fastq_LM'
	# use_strain = None

	print('Creating summaries')

	os.chdir(baseDir)

	run_data = {}
	for s in glob.glob(pjoin(sequenceDir,'*_R1_*.fastq.gz')):
		sample = SequencingSample()
		sample.get_DataFromReadPairs(read1_filename=s)
		sample.check_ReadExistence(directory=sequenceDir)
		sample.get_SampleDirPath(directory=baseDir)

		print(sample.samplename)

		sample_data = {}
		# check errors in log
		if os.path.exists(pjoin(baseDir,'sampleLogs',sample.samplename+'.txt')):
			time_ran = None
			with open(pjoin(baseDir,'sampleLogs',sample.samplename+'.txt'), 'r') as infile:
				for l in infile:
					if l.find('FluPipeline Error') > -1:
						sample_data['errors'] = l.rstrip('\n')
					else:
						sample_data['errors'] = None
					if l.find('FluPipeline Time:') > -1:
						time_ran = l[l.find('me: ')+4:].replace('\n','') #use 'me: ' as an anchor + 4 positions to extract just the time
		if time_ran == None:
			sample_data['run_time'] = None
		else:
			sample_data['run_time'] = time_ran.split('.')[0]+' H:M:S'

		sample_data['reference_selection'] = {}
		## reference strain assigned ## ---------------
		if use_strain == None:
			try:
				reference_strain = select_BestReference(samplename=sample.samplename, directory=sample.dirpath)
			except:
				reference_strain = None
			sample_data['reference_selection']['best_reference_strain'] = reference_strain
		else:
			reference_strain = use_strain
			sample_data['reference_selection']['best_reference_strain'] = reference_strain


		## reference_coverage is a dataframe ## ---------------
		if use_strain == None:
			df = safe_ReadCSV(pjoin(sample.dirpath,'reference_coverage_'+sample.samplename+'.csv'))
			if len(df) == 0:
				sample_data['reference_selection']['coverage_statistics'] = None
			else:
				df = df[['reference','segment_count','segment_coverage_average','segment_average_read_depth']]
				df = df.rename(columns={'segment_coverage_average':'avg_coverage','segment_average_read_depth':'avg_depth'})
				df['reference'] = [i.replace('_coverage_stats.csv','') for i in df.reference.tolist()]
				sample_data['reference_selection']['coverage_statistics'] = df.to_dict()

		##  read_stats is a dataframe ## ---------------
		df = safe_ReadCSV(pjoin(sample.dirpath,'fastp_stats_'+sample.samplename+'.csv'))
		if len(df) == 0:
			sample_data['read_statistics'] = None
		else:
			df = df.drop(columns='sample')
			sample_data['read_statistics'] = df.to_dict()

		## coverage_stats is a dataframe ## ---------------
		df = safe_ReadCSV(pjoin(sample.dirpath,sample.samplename+'_covstats.txt'), sep_='\t')
		if len(df) == 0:
			sample_data['coverage_statistics'] = None
		else:
			df = df.rename(columns={df.columns[0]:'ID'})
			df = df[['ID','Avg_fold','Median_fold','Length','Covered_bases']]
			df['Covered_proportion'] = (df['Covered_bases']/df['Length']).round(2)
			df = df.drop(columns=['Covered_bases'])
			sample_data['coverage_statistics'] = df.to_dict()

		## all_variants is a dataframe ## ---------------
		df = safe_ReadCSV(pjoin(sample.dirpath,sample.samplename+'_allVariants.csv'))
		if len(df) == 0:
			sample_data['all_variants'] = None
		else:
			df = format_VariantTable(df)
			sample_data['all_variants'] = df.to_dict()

		## passing_variants is a dataframe ## ---------------
		df = safe_ReadCSV(pjoin(sample.dirpath,sample.samplename+'_passingVariants.csv'))
		if len(df) == 0:
			sample_data['passing_variants'] = None
		else:
			df = format_VariantTable(df)
			sample_data['passing_variants'] = df.to_dict()

		## major_variants is a dataframe ## ---------------
		df = safe_ReadCSV(pjoin(sample.dirpath,sample.samplename+'_majorVariants.csv'))
		if len(df) == 0:
			sample_data['major_variants'] = None
		else:
			df = format_VariantTable(df)
			sample_data['major_variants'] = df.to_dict()

		##### IVAR ####

		sample_data['ivar'] = {}

		## all_variants is a dataframe ## ---------------
		df = safe_ReadCSV(pjoin(sample.dirpath,sample.samplename+'_ivar',sample.samplename+'_ivar_allVariants.csv'))
		if len(df) == 0:
			sample_data['ivar']['all_variants'] = None
		else:
			df = format_VariantTable(df)
			sample_data['ivar']['all_variants'] = df.to_dict()

		## passing_variants is a dataframe ## ---------------
		df = safe_ReadCSV(pjoin(sample.dirpath,sample.samplename+'_ivar',sample.samplename+'_ivar_passingVariants.csv'))
		if len(df) == 0:
			sample_data['ivar']['passing_variants'] = None
		else:
			df = format_VariantTable(df)
			sample_data['ivar']['passing_variants'] = df.to_dict()

		## major_variants is a dataframe ## ---------------
		df = safe_ReadCSV(pjoin(sample.dirpath,sample.samplename+'_ivar',sample.samplename+'_ivar_majorVariants.csv'))
		if len(df) == 0:
			sample_data['ivar']['major_variants'] = None
		else:
			df = format_VariantTable(df)
			sample_data['ivar']['major_variants'] = df.to_dict()

		# save json report
		with open(pjoin(sample.dirpath,sample.samplename+'_report.json'), 'w') as fp:
			json.dump(sample_data, fp, indent=4)


		# aggregate run-wide data 
		run_data[sample.samplename] = {}
		run_data[sample.samplename]['errors'] = sample_data['errors']
		if sample_data['run_time'] == None:
			run_data[sample.samplename]['run_time'] = sample_data['run_time']
		else:
			run_data[sample.samplename]['run_time'] = sample_data['run_time'].replace(' H:M:S','')
		run_data[sample.samplename]['best_reference_strain'] = sample_data['reference_selection']['best_reference_strain']
		run_data[sample.samplename]['reads_after_filtering'] = pd.DataFrame(sample_data['read_statistics'])['total_reads_after'].item()

		df = pd.DataFrame(sample_data['coverage_statistics'])

		if len(df) > 0:
			run_data[sample.samplename]['segments_covered'] = len(df)
			run_data[sample.samplename]['segments_covered_90'] = len(df[df['Covered_proportion'] >= .9])
			run_data[sample.samplename]['avg_prop_coverage'] = round(df['Covered_proportion'].mean(),2)
			run_data[sample.samplename]['avg_fold_coverage'] = round(df['Avg_fold'].mean(),2)
		else:
			run_data[sample.samplename]['segments_covered'] = int(0)
			run_data[sample.samplename]['segments_covered_90'] = float(0)
			run_data[sample.samplename]['avg_prop_coverage'] = float(0)
			run_data[sample.samplename]['avg_fold_coverage'] = float(0)

		del df

	df = pd.DataFrame(run_data).T.sort_values(['best_reference_strain','errors'])

	df.to_csv(pjoin(baseDir,'runReport.csv'),index=True)

	with open(pjoin(baseDir,'runReport.json'), 'w') as fp:
		json.dump(run_data, fp, indent=4)








