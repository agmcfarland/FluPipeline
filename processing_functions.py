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
from processing_classes import SequencingSample, RunLogger, SampleTracker

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


def flu_Pipeline(
	baseDirectory,
	Rscript,
	softwareDir,
	sample,
	pipeline_used,
	sequenceDataDir,
	referenceStrainsDir,
	force_overwrite,
	BWA_path,
	samtoolsbin_path,
	bcftoolsbin_path):
	'''
	Runs through all steps of the flu pipeline. 
	'''

	# ## troubleshooting inputs start ##
	# baseDirectory='/home/agmcfarland/flu_project/test/test4'
	# Rscript='Rscript'
	# softwareDir = '/home/agmcfarland/flu_project/flu_pipeline'
	# sample = SequencingSample()
	# sample.get_DataFromReadPairs(read1_filename='H1N1pdm_ref_snp_R1_001.fastq.gz')
	# sample.get_SampleDirPath(directory=baseDirectory)
	# pipeline_used='snp'
	# sequenceDataDir = '/home/agmcfarland/flu_project/test/test_data'
	# referenceStrainsDir = '/home/agmcfarland/flu_project/test/test_data' #directory reference strains to find the best match for a given readset
	# pipeline_used = 'snp' #options: snp, phylo
	# process_method = 'bushman_artic_v2'
	# force_overwrite = True #deletes and remakes the sample folder in sampleOutputs
	# BWA_path =  '/home/agmcfarland/miniconda3/envs/FluPipeline_env/bin/bwa'#'/home/everett/ext/bwa' #to keep bwa version consistent
	# samtoolsbin_path = '/home/agmcfarland/miniconda3/envs/FluPipeline_env/bin'#'/home/everett/ext/samtools/bin' #to keep samtools verison consistent
	# bcftoolsbin_path = '/home/agmcfarland/miniconda3/envs/FluPipeline_env/bin'
	# ## troubleshooting inputs end ##


	## ==================================Prepare sample run directory================================================
	## ==============================================================================================================
	# set and create directory to store all intermediate data in sample-specific folders. Switch to that directory afterwards.

	if force_overwrite==True:
		try: 
			shutil.rmtree(sample.dirpath)
		except:
			pass

	os.makedirs(sample.dirpath, exist_ok=True)
	os.chdir(sample.dirpath)

	start_run_timer = datetime.now()
	sample_logger = RunLogger(directory=pjoin(baseDirectory,'sampleLogs'),filename=sample.samplename)
	sample_logger.delete_LogFile()
	sample_logger.initialize_Log()
	sample_logger.add_Message('Starting run')
	sample_logger.add_Message('sample: {}'.format(sample.samplename))

	## ==================================Find and assign reference strain============================================
	## ==============================================================================================================
	try:
		sample_logger.add_Message('Assigning reference strain')
		# copy reference fastas to sample directory 
		reference_fasta_list = glob.glob(pjoin(referenceStrainsDir,'*.fasta')) #glob of fasta files
		[shutil.copy(rf,sample.dirpath) for rf in reference_fasta_list] # copy to sample dir


		# quality, adaptor trimming with fastp
		subprocess.run([
			'fastp', 
			'--in1', pjoin(sequenceDataDir,sample.read1_filename), #readfileR1
			'--in2', pjoin(sequenceDataDir,sample.read2_filename), #readfileR2
			'--out1', 'fastp_trimmed_{}'.format(sample.read1_filename), 
			'--out2', 'fastp_trimmed_{}'.format(sample.read2_filename),
			'-j', 'fastp_stats_{}'.format(sample.samplename), #json output
			'-q', '30'#, #quality score 30
			#'--disable_adapter_trimming' # adaptor trimming when adaptors are already trimmed can lead to unwanted errors
			])

		# randomly select 1000 reads with reformat.sh
		subprocess.run([
			'reformat.sh',
			'in1=fastp_trimmed_{}'.format(sample.read1_filename), 
			'in2=fastp_trimmed_{}'.format(sample.read2_filename),
			'out1=subset_fastp_trimmed_{}'.format(sample.read1_filename), 
			'out2=subset_fastp_trimmed_{}'.format(sample.read2_filename),
			'samplereadstarget=2000'
		])


		# for each reference available get coverage stats
		for reference in glob.glob('*.fasta'):
			calculate_ReferenceCoverage(sequenceDataDir=sequenceDataDir,reference=reference, BWA_path=BWA_path, samtoolsbin_path=samtoolsbin_path, read1_filename=sample.read1_filename, read2_filename=sample.read2_filename,samplename=sample.samplename)

		# select the best reference strain and assign it to variables used by assemble.R
		summarize_ReferenceCoverage(samplename=sample.samplename)
		summarize_HACoverage(samplename=sample.samplename)
		reference_strain = select_BestReference(samplename=sample.samplename).replace('.fasta_coverage_stats.csv','.fasta')
		refGenomeFasta = reference_strain
		refGenomeBWA = reference_strain

	except:
		sample_logger.add_Message('failure at strain assignment',level='warning')
		end_run_timer = datetime.now()-start_run_timer
		return(pd.DataFrame({'sample':[sample.samplename], 'processed':[pipeline_used], 'date_run':[datetime.now()], 'runtime':[end_run_timer.seconds], 'info':['failure at strain assignment']}))


	## ==================================assemble.R find variants====================================================
	## ==============================================================================================================
	try:
		sample_logger.add_Message('Finding variants')
		# run assemble.R 
		os.system('{} {}/assemble.R --outputFile {} --workDir {} --process_method {} --softwareDir {} --R1 {} --R2 {} --refGenomeFasta {} --refGenomeBWA {} --bwaPath {} --samtoolsBin {} --bcftoolsBin {}'.
		format(
			Rscript, 
			softwareDir,
			pjoin(sample.dirpath,sample.samplename+'.Rdata'), #--outputFile
			pjoin(sample.dirpath), #--workDir
			'bushman_artic_v2', #--process_method
			softwareDir, #--softwareDir
			pjoin(sequenceDataDir,sample.read1_filename), #--R1
			pjoin(sequenceDataDir,sample.read2_filename), #--R2
			refGenomeFasta, #--refGenomeFasta
			refGenomeBWA, #--refGenomeBWA
			BWA_path,
			samtoolsbin_path,
			bcftoolsbin_path
		))

	except:
		sample_logger.add_Message('failure at variant calling',level='warning')
		end_run_timer = datetime.now()-start_run_timer
		return(pd.DataFrame({'sample':[sample.samplename], 'processed':[pipeline_used], 'date_run':[datetime.now()], 'runtime':[end_run_timer.seconds], 'info':['failure at variant calling']}))

	try:
		sample_logger.add_Message('Generating sample report')
		# generate sample report
		os.system('{} {}/report_runner.R --softwareDir {} --report_type {} --sampleDir {} --samplename {}'.
			format(
				Rscript,
				softwareDir,
				softwareDir, #--softwareDir
				'sample', #--report_type
				sample.dirpath, #--sampleDir
				sample.samplename #--samplename
				))

	except:
		sample_logger.add_Message('failure at sample report generation',level='warning')
		end_run_timer = datetime.now()-start_run_timer
		return(pd.DataFrame({'sample':[sample.samplename], 'processed':[pipeline_used], 'date_run':[datetime.now()], 'runtime':[end_run_timer.seconds], 'info':['failure at sample report generation']}))

	# ==================================clean up and record run=====================================================
	# ==============================================================================================================
	# sample.cleanup_OutputFiles(level='intermediate')
	sample_logger.add_Message('Finished processing sample')
	end_run_timer = datetime.now()-start_run_timer
	# record the end of the run 
	return(pd.DataFrame({'sample':[sample.samplename], 'processed':[pipeline_used], 'date_run':[datetime.now()], 'runtime':[end_run_timer.seconds], 'info':['none']}))


def calculate_ReferenceCoverage(sequenceDataDir, reference, BWA_path, samtoolsbin_path, read1_filename, read2_filename, samplename):
	'''
	Runs through all steps of a read mapping pipeline and calculates the coverage and depth for a given read file reference combo.
	'''

# ## troubleshooting inputs start ##
# read1_filename = sample.read1_filename
# read2_filename = sample.read2_filename
# samplename = sample.samplename
# reference = 'H1N1pdm_ref.fasta'
# ## troubleshooting inputs end ##

	print('ref: ',reference)
	print('read1_filename: ',read1_filename)
	print('read2_filename: ',read2_filename)
	print('samplename: ',samplename)

	## index reference
	subprocess.run([
		BWA_path,
		'index',
		reference
		])

	## align random read subset to reference with bwa
	os.system('{} mem -M {} subset_fastp_trimmed_{} subset_fastp_trimmed_{} > {}_bwa_alignment_{}.sam'.format(BWA_path, reference, read1_filename, read2_filename, reference, samplename))

	## sam to bam and sort and index alignment using samtools
	os.system('{}/samtools view -S -b {}_bwa_alignment_{}.sam > {}_bwa_alignment_{}.bam'.format(samtoolsbin_path,reference,samplename,reference,samplename))
	os.system('{}/samtools sort -o {}_bwa_alignment_{}.bam {}_bwa_alignment_{}.bam'.format(samtoolsbin_path,reference,samplename,reference,samplename))
	os.system('{}/samtools index {}_bwa_alignment_{}.bam'.format(samtoolsbin_path,reference,samplename))

	## pileup data with samtools
	os.system('{}/samtools mpileup -A -a -Q 0 -o {}_pileup_{}.txt -d 10000 -f {} {}_bwa_alignment_{}.bam'.format(samtoolsbin_path, reference,samplename, reference, reference,samplename))

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
	df_cov.to_csv(pjoin('reference_coverage_{}.csv'.format(samplename)), index=False)

	print(df_cov)


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
	df_ha.to_csv(pjoin('reference_ha_coverage_{}.csv'.format(samplename)), index=False)



def select_BestReference(samplename):
	'''
	Selects the best reference from the output of summarize_ReferenceCoverage
	(most segments, highest coverage, highest read depth). 
	'''
	df = pd.read_csv(pjoin('reference_ha_coverage_{}.csv'.format(samplename)))
	df = df.sort_values(['ha_coverage','ha_read_depth'], ascending = [False,False])
	return(df.head(1)['reference'].item())



