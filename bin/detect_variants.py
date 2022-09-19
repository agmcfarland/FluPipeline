
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

def write_ConsensusSeq(filename, genome_consensus):
	'''
	'''
	with open(filename, 'w') as outfile:
		for chrom, seq in genome_consensus.items():
			outfile.write(f'>{chrom}\n')
			[outfile.write(i) if e%60!=0 else outfile.write(f'{i}\n') for e, i in enumerate(seq, 1)]
			outfile.write('\n')

def is_NaN(x):
	'''
	'''
	return x != x

def check_Vcf(df_vcf):
	'''
	'''
	if len(df_vcf) == 0:
		raise Warning('No variants found')
		remove_IntermediateFiles()
		sys.exit()


def remove_IntermediateFiles():
	'''
	'''
	### Remove files
	pass
	# remove bam files
	for f in glob.glob('*.bam'):
		if f == f'{samplename}_genome.indels.filt.qual.sorted.bam':
			continue
		elif f == f'{samplename}_genome.indels.filt.qual.sorted.bam.bai':
			continue
		elif f == f'{samplename}_genome.filt.qual.sorted.bam':
			continue
		elif f == f'{samplename}_genome.filt.qual.sorted.bam.bai':
			continue
		os.remove(f)
	# remove index files
	for f in os.listdir():
		if f.endswith('.sa'):
			os.remove(f)
		elif f.endswith('.pac'):
			os.remove(f)
		elif f.endswith('.bwt'):
			os.remove(f)
		elif f.endswith('.ann'):
			os.remove(f)
		elif f.endswith('.amb'):
			os.remove(f)
		else:
			continue


if __name__ == '__main__':
	args = sys.argv[1:]
	parser = argparse.ArgumentParser(prog='detect_variants')
	parser.add_argument('--baseDir', type=str, default=os.getcwd(), help='none', metavar='none')
	parser.add_argument('--logDir', type=str, default=os.getcwd(), help='none', metavar='none')
	parser.add_argument('--variant_caller', type=str, default=os.getcwd(), help='none', metavar='none')
	parser.add_argument('--R1', type=str, default=None, help='none', metavar='none')
	parser.add_argument('--R2', type=str, default=None, help='none', metavar='none')
	parser.add_argument('--refGenomeFasta', type=str, default=None, help='none', metavar='none')
	parser.add_argument('--minAmpliconLength', type=int, default=50, help='none', metavar='none')
	parser.add_argument('--maxAmpliconLength', type=int, default=350, help='none', metavar='none')
	parser.add_argument('--removeNTsFromAlignmentEnds', type=int, default=3, help='none', metavar='none')
	parser.add_argument('--BWAmappingScore', type=int, default=60, help='none', metavar='none')
	parser.add_argument('--minVariantPhredScore', type=int, default=20, help='none', metavar='none')
	parser.add_argument('--majorIndelVariantThreshold', type=float, default=0.8, help='none', metavar='none')
	parser.add_argument('--majorVariantThreshold', type=float, default=0.5, help='none', metavar='none')
	parser.add_argument('--minorVariantThreshold', type=float, default=0.05, help='none', metavar='none')
	parser.add_argument('--consensus_sequence', action='store_true', default=False, help='none')
	parser.add_argument('--consensus_masking_threshold', type=int, default=1, help='none', metavar='none')
	parser.add_argument('--minimum_read_depth', type=int, default=10, help='none', metavar='none')

	## ADDED ##
	# baseDir = '/data/flu_project/benchmarking_project/output'
	# variant_caller = 'bcftools'

	# ## OG ## 
	## Troubleshooting inputs
	samplename = 'CHOA-063_S38'
	baseDir = '/data/flu_project/benchmarking_project/compare_caller/lofreq/sampleOutputs/CHOA-063_S38'
	logDir = baseDir
	R1 = 'fastp_trimmed_CHOA-063_S38_R1_001.fastq.gz'
	R2 = 'fastp_trimmed_CHOA-063_S38_R2_001.fastq.gz'
	refGenomeFasta = 'H3N2_ref.fasta'
	minAmpliconLength = 50
	maxAmpliconLength = 350
	minVariantPhredScore = 30
	removeNTsFromAlignmentEnds = 3
	BWAmappingScore = 60
	minorVariantThreshold = 0.05
	majorVariantThreshold = 0.8
	majorIndelVariantThreshold = 0.8
	variant_caller = 'bcftools'
	minimum_read_depth = 10

	args = parser.parse_args()

	baseDir=args.baseDir
	logDir=args.logDir
	variant_caller=args.variant_caller
	R1=args.R1
	R2=args.R2
	refGenomeFasta=args.refGenomeFasta
	minAmpliconLength=args.minAmpliconLength
	maxAmpliconLength=args.maxAmpliconLength
	removeNTsFromAlignmentEnds=args.removeNTsFromAlignmentEnds
	BWAmappingScore=args.BWAmappingScore
	minVariantPhredScore=args.minVariantPhredScore
	majorIndelVariantThreshold=args.majorIndelVariantThreshold
	majorVariantThreshold=args.majorVariantThreshold
	minorVariantThreshold=args.minorVariantThreshold
	consensus_sequence=args.consensus_sequence
	consensus_masking_threshold=args.consensus_masking_threshold
	minimum_read_depth=args.minimum_read_depth


	if os.path.exists(R1) == False:
		raise ValueError(f'R1 does not exist:\n{R1}')
	if os.path.exists(R2) == False:
		raise ValueError(f'R2 does not exist:\n{R2}')
	if os.path.exists(refGenomeFasta) == False:
		raise ValueError(f'Reference genome does not exist:\n{refGenomeFasta}')

	# define sample name and specify explicit path
	os.chdir(baseDir)
	samplename = os.path.basename(os.getcwd())

	# start logging
	logger = RunLogger(directory=logDir,filename=samplename)
	logger.initialize_FileHandler()
	logger.logger.info('Running detect_variants.py\n')
	logger.logger.info('Processing sample: {}'.format(samplename))

	
	logger.logger.info('Program arguments:\n')
	arguments_list = vars(args)
	for k,v in arguments_list.items():
		logger.logger.info('{}: {}\n'.format(k,v))

	### Prep BAM ####

	# copy reference genome fasta to the working directory to use for indexing
	logger.logger.info('Making BWA index...')
	if os.path.exists(pjoin(os.getcwd(),os.path.basename(refGenomeFasta))) == False:
		shutil.copy(refGenomeFasta, os.getcwd())

	refGenomeFasta = os.path.basename(refGenomeFasta)
	# make index
	call_Command(cmd=
		[
		'bwa','index',refGenomeFasta
		],
		logger_=logger)

	# align reads to reference. make a sam file
	logger.logger.info('Aligning reads...')
	call_Command(cmd=
		f'bwa mem -M {refGenomeFasta} {R1} {R2} > {samplename}_genome.sam'
		,
		logger_=logger,
		shell_=True)

	# convert sam to bam and delete sam file (save space!)
	call_Command(cmd=
		f'samtools view -S -b -h {samplename}_genome.sam > {samplename}_genome.bam'
		,
		logger_=logger,
		shell_=True)
	os.remove(f'{samplename}_genome.sam')

	# Remove non-paired mates and filter alignment by insert size
	call_Command(cmd=
		f"samtools view -F 8 -h {samplename}_genome.bam | awk 'length($10) > {minAmpliconLength} && length($10) < {maxAmpliconLength} || $1 ~ /^@/' | samtools view -b -h > {samplename}_genome.filt.bam"
		,
		logger_=logger,
		shell_=True)

	# trim ends of alignments
	# use trimBam command and then rename the output file to be suffixed with _genome.filt.bam
	if removeNTsFromAlignmentEnds > 0:
		call_Command(cmd=
			f'bam trimBam {samplename}_genome.filt.bam {samplename}_genome.filt.endTrim.bam {removeNTsFromAlignmentEnds} --clip'
			,
			logger_=logger,
			shell_=True)
		os.rename(f'{samplename}_genome.filt.endTrim.bam', f'{samplename}_genome.filt.bam')

	# Remove read pairs with mapping qualities below the provided min score
	call_Command(cmd=
		f'samtools view -h -q {BWAmappingScore} -b {samplename}_genome.filt.bam > {samplename}_genome.filt.qual.bam'
		,
		logger_=logger,
		shell_=True)

	# Get number of aligned reads and check that it is greater than 0. Record error. 
	lengthAlignedReadsIDs = len(subprocess.Popen(f'samtools view {samplename}_genome.filt.qual.bam | cut -f 1 | uniq', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0].decode('utf-8'))
	if lengthAlignedReadsIDs == 0:
		logger.logger.exception('FluPipeline Error: No reads mapped after filtering.\n')
		sys.exit()
	logger.logger.info(f'{lengthAlignedReadsIDs} reads aligned\n')

	# Sort and index filtered genome alignment.
	logger.logger.info(f'Sorting and and indexing genome alignment...\n')
	call_Command(cmd=
		f'samtools sort -o {samplename}_genome.filt.qual.sorted.bam {samplename}_genome.filt.qual.bam'
		,
		logger_=logger,
		shell_=True)
	call_Command(cmd=
		f'samtools index {samplename}_genome.filt.qual.sorted.bam'
		,
		logger_=logger,
		shell_=True)


	# Pileup reads. exit if reads failed to pileup.
	logger.logger.info('Piling up reads to get coverage...\n')
	call_Command(cmd=
		f'pileup.sh in={samplename}_genome.filt.qual.sorted.bam ref={refGenomeFasta} out={samplename}_covstats.txt basecov={samplename}_basecov.txt countgc=f overwrite=t 32bit=t minmapq={BWAmappingScore}'
		,
		logger_=logger,
		shell_=True)
	if len(pd.read_csv(f'{samplename}_covstats.txt', sep='\t')) == 0:
		logger.logger.exception('FluPipeline Error: No reads pileup..\n')


	### Call variants ####

	logger.logger.info(f'Calling variants with {variant_caller}...\n')

	if variant_caller == 'bbtools':
		call_Command(cmd=
		f'callvariants.sh in={samplename}_genome.filt.qual.sorted.bam ref={refGenomeFasta} out={samplename}_allVariants.vcf shist={samplename}_variantQualityHisto.txt rarity=0 overwrite=t' #clearfilters
		,
		logger_=logger,
		shell_=True)

		df_vcf = pd.read_csv(f'{samplename}_allVariants.vcf', sep='\t',comment='#',header=None,
			names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','OTHER'])

		check_Vcf(df_vcf = df_vcf)


		vcf_info = df_vcf['INFO'].str.split(';',expand=True)

		vcf_info_cols = [i[:i.find('=')] for i in vcf_info.iloc[0,:].tolist()]

		vcf_info[vcf_info.columns] = vcf_info.apply(lambda x: [i[i.find('=')+1:] for i in x])

		vcf_info.columns = vcf_info_cols

		# Rename to TYPE for consistency across variant callers
		vcf_info.rename(columns={'TYP':'TYPE'}, inplace=True)

		df_vcf = pd.concat([df_vcf, vcf_info], axis=1)

		df_vcf['AF'] = df_vcf['AF'].astype(float)

		df_vcf['QUAL'] = df_vcf['QUAL'].astype(float)

		# write to file
		df_vcf.to_csv(f'{samplename}_allVariants.csv', index=None)

		# passing variants
		df_vcf = df_vcf[(df_vcf['AF']>=minorVariantThreshold) & (df_vcf['QUAL']>=float(minVariantPhredScore))]
		df_vcf.to_csv(f'{samplename}_passingVariants.csv', index=None)

		# major variants
		df_vcf = df_vcf[df_vcf['AF']>0.5]
		df
		df_vcf.to_csv(f'{samplename}_majorVariants.csv', index=None)


	if variant_caller == 'lofreq':
		# add indel qualities to bam. save to new filename
		call_Command(cmd=
		f'lofreq indelqual --dindel -f {refGenomeFasta} --out {samplename}_genome.indels.filt.qual.sorted.bam {samplename}_genome.filt.qual.sorted.bam'
		,
		logger_=logger,
		shell_=True)

		# index indel bam file
		call_Command(cmd=
		f'samtools index {samplename}_genome.indels.filt.qual.sorted.bam'
		,
		logger_=logger,
		shell_=True)

		# call snps 
		call_Command(cmd=
		f'lofreq call-parallel --pp-threads 2  -N --call-indels -f {refGenomeFasta} -o {samplename}_lofreq.vcf {samplename}_genome.indels.filt.qual.sorted.bam'
		,
		logger_=logger,
		shell_=True)

		# Rework VCF into table with columns for specific INFO sections
		df_vcf = pd.read_csv(f'{samplename}_lofreq.vcf', sep='\t',comment='#',header=None,
			names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'])

		check_Vcf(df_vcf = df_vcf)

		vcf_info = df_vcf['INFO'].str.split(';',expand=True)

		vcf_info.rename(columns={
			vcf_info.columns[0]: 'DP',
			vcf_info.columns[1]: 'AF',
			vcf_info.columns[2]: 'SB',
			vcf_info.columns[3]: 'DP4'}, inplace = True)

		if len(vcf_info.columns) == 6:
			vcf_info.rename(columns={
			vcf_info.columns[4]: 'TYPE',
			vcf_info.columns[5]: 'INDEL_LENGTH'}, inplace = True)
			vcf_info['TYPE'] = ['SUB' if i is None else i for i in vcf_info.TYPE.tolist()]
			vcf_info['INDEL_LENGTH'] = vcf_info['INDEL_LENGTH'].str.replace('HRUN=','')
			vcf_info['INDEL_LENGTH'] = [0 if i is None else i for i in vcf_info.INDEL_LENGTH.tolist()]
			vcf_info['INDEL_LENGTH'] = vcf_info['INDEL_LENGTH'].astype(int)

		else:
			# add TYPE column
			vcf_info['TYPE'] == 'SUB'
			vcf_info['INDEL_LENGTH'] = 0

		vcf_info['DP'] =  vcf_info['DP'].str.replace('DP=','').astype(int)
		vcf_info['AF'] =  vcf_info['AF'].str.replace('AF=','').astype(float)
		vcf_info['SB'] =  vcf_info['SB'].str.replace('SB=','').astype(int)
		vcf_info['DP4'] =  vcf_info['DP4'].str.replace('DP4=','')

		df_vcf = pd.concat([df_vcf, vcf_info], axis=1)

		# write to file
		df_vcf.to_csv(f'{samplename}_allVariants.csv', index=None)

		df_vcf = df_vcf[(df_vcf['AF']>=minorVariantThreshold) & (df_vcf['QUAL']>=float(minVariantPhredScore)) & (df_vcf['DP']>minimum_read_depth)]
		df_vcf.to_csv(f'{samplename}_passingVariants.csv', index=None)

		df_vcf = df_vcf[df_vcf['AF']>0.5]
		df_vcf.to_csv(f'{samplename}_majorVariants.csv', index=None)


		# get segment names and sizes to add to vcf
		segment_name_size = {}
		for record in SeqIO.parse(refGenomeFasta,'fasta'):
			segment_name_size[record.id] = len(str(record.seq))

		# make vcf from filter
		with open(f'{samplename}_lofreq.vcf', 'r') as infile:
			with open(f'consensus_input.vcf','w') as outfile:
				header_line = True
				for l in infile:
					if header_line == True:
						if l.find('#CHROM') == -1:
							outfile.write(l)
							if l.find('##reference=') > -1:
								for k,v in segment_name_size.items():
									outfile.write(f'##contig=<ID={k},length={v}>\n')
							if l.find('Homopolymer') > -1:
								outfile.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n') #['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'] #https://github.com/CSB5/lofreq/blob/master/src/tools/scripts/lofreq2_add_fake_gt.py
						else:
							l = l.replace('\n','\tFORMAT\n')
							outfile.write(l)
							header_line = False


		vcf_filtered = pd.read_csv(f'{samplename}_majorVariants.csv')[['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO']]
		vcf_filtered['FORMAT'] = 'GT'
		vcf_filtered.to_csv('consensus_vcf_table.txt',sep='\t',header=None,index=None)
		os.system(f'cat consensus_vcf_table.txt >> consensus_input.vcf')

		os.system('head -n 50 consensus_input.vcf')
		os.system(f'bcftools norm -f {refGenomeFasta} -o consensus_input_norm.vcf consensus_input.vcf')
		os.system(f'bgzip -c consensus_input.vcf > consensus_input.vcf.gz')
		os.system(f'tabix consensus_input.vcf.gz')
		os.system(f'bcftools consensus --mark-ins lc --mark-snv lc --mark-del lc -f {refGenomeFasta} -o bcftools_consensus.fasta consensus_input.vcf.gz')

	if variant_caller == 'bcftools':
		# -A: count orphans
		# -a: output all positions, including those with 0 depth
		# -Q: minimum base quality for a base to be considered 
		# -d max number of reads to considerper base
		# -m: multiallelic-caller
		# -v variant sites only
		# -p: variant if pvalue greater than
		call_Command(cmd=
		f'bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR -d 1000000 -Q 0 -f {refGenomeFasta} {samplename}_genome.filt.qual.sorted.bam | bcftools call -p 0.01 --ploidy 1 -m -v -Ov -o {samplename}_bcftools.vcf'
		,
		logger_=logger,
		shell_=True)

		df_vcf = pd.read_csv(f'{samplename}_bcftools.vcf', sep='\t',comment='#',header=None,
			names=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO', 'FORMAT', 'OTHER'])

		info_dict = {}
		for e, row in df_vcf.iterrows():
			info = ';'+row.INFO
			info_split  = info.split(';')[1:]
			temp_dict = {}
			TYPE = 'SUB'
			for i in info_split:
				if i == 'INDEL':
					TYPE = 'INDEL'
					continue
				else:
					temp_dict[i.split('=')[0]] = i.split('=')[1] 
			temp_dict['TYPE'] = TYPE
			info_dict[e] = temp_dict

		vcf_info = pd.DataFrame(info_dict).T

		df_vcf = pd.concat([df_vcf, vcf_info], axis=1)

		check_Vcf(df_vcf = df_vcf)

		df_vcf['AD_ORIG'] = df_vcf['AD']
		df_vcf['AD'] = [i.split(',')[1] for i in df_vcf['AD'].tolist()]
		df_vcf['AD'] = df_vcf['AD'].astype(float)
		df_vcf['DP'] = df_vcf['DP'].astype(float)
		df_vcf['MQ'] = df_vcf['MQ'].astype(int)
		df_vcf['AF'] = df_vcf['AD']/df_vcf['DP']

		df_vcf['AF'] = df_vcf['AF'].astype(float)

		df_vcf['QUAL'] = df_vcf['QUAL'].astype(float)

		# write to file
		df_vcf.to_csv(f'{samplename}_allVariants.csv', index=None)

		# passing variants
		df_vcf = df_vcf[(df_vcf['AF']>=minorVariantThreshold) & (df_vcf['QUAL']>=float(minVariantPhredScore)) & (df_vcf['DP']>minimum_read_depth) & (df_vcf['MQ']>=BWAmappingScore)]
		df_vcf.to_csv(f'{samplename}_passingVariants.csv', index=None)

		# major variants
		df_vcf = df_vcf[((df_vcf['TYPE']=='SUB')&(df_vcf['AF']>majorVariantThreshold))| ((df_vcf['TYPE']=='INDEL') & (df_vcf['AF']>majorIndelVariantThreshold))]

		df_vcf.to_csv(f'{samplename}_majorVariants.csv', index=None)


		# make vcf from filter
		with open(f'{samplename}_bcftools.vcf', 'r') as infile:
			with open(f'consensus_input.vcf','w') as outfile:
				header_line = True
				for l in infile:
					if header_line == True:
						outfile.write(l)
						if l.find('#CHROM') > -1:
							header_line = False
		vcf_filtered = pd.read_csv(f'{samplename}_majorVariants.csv')[['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO', 'FORMAT', 'OTHER']]
		vcf_filtered.to_csv('consensus_vcf_table.txt',sep='\t',header=None,index=None)
		os.system(f'cat consensus_vcf_table.txt >> consensus_input.vcf')
		os.system(f'bcftools norm -f {refGenomeFasta} -o consensus_input_norm.vcf consensus_input.vcf')
		os.system(f'bgzip -c consensus_input.vcf > consensus_input.vcf.gz')
		os.system(f'tabix consensus_input.vcf.gz')
		os.system(f'bcftools consensus --mark-ins lc --mark-snv lc --mark-del lc -f {refGenomeFasta} -o bcftools_consensus.fasta consensus_input.vcf.gz')



	### Make consensus sequence ####

	if consensus_sequence == True:
		logger.logger.info('Writing consensus sequence...\n')
		# make consensus sequence
		df_vcf = pd.read_csv(f'{samplename}_majorVariants.csv')

		df_reference = pd.DataFrame()
		for record in SeqIO.parse(refGenomeFasta,'fasta'):
			temp = pd.DataFrame([[i, e+1] for e, i in enumerate(str(record.seq))])
			temp['CHROM'] = record.id
			if len(df_reference) == 0:
				df_reference = temp
			else:
				df_reference = pd.concat([df_reference, temp], ignore_index=True)

		df_reference.columns = ['REF','POS', 'CHROM']
		df_reference = df_reference[['CHROM','POS','REF']]

		df_pileup = pd.read_csv(f'{samplename}_basecov.txt', sep='\t', names=['CHROM','POS','PILEUP_COV'], comment='#')
		df_pileup['POS'] = df_pileup['POS']+1

		df_reference = df_reference.merge(df_pileup, on=['CHROM', 'POS'])

		df_reference = df_reference.merge(df_vcf[['CHROM','POS','REF','ALT','TYPE']], on=['CHROM','POS'], how='left')

		chrom_order =  [i.split('_')[2] for i in df_reference['CHROM'].unique().tolist()]

		df_reference['chrom_order'] = df_reference['CHROM'].apply(lambda x: x.split('_')[2])

		df_reference = df_reference.sort_values(['chrom_order','POS'], ascending=[True, True])

		chromosome_list = [i for i in df_reference['CHROM'].unique().tolist()]

		genome_consensus = {'unmasked':{}, 'masked':{}}
		for chromosome in chromosome_list:
			consensus_seq = ''
			masked_consensus_seq = ''
			for e, row in df_reference[df_reference['CHROM']==chromosome].iterrows():
				base = ''
				if is_NaN(row.TYPE) == False:
					if row.TYPE == 'SUB':
						base = row.ALT
				else:
					base = row.REF_x
				consensus_seq += base

				if row.PILEUP_COV < minimum_read_depth:
					base = 'N'
				masked_consensus_seq += base

			genome_consensus['unmasked'][chromosome] = consensus_seq
			genome_consensus['masked'][chromosome] = masked_consensus_seq

		write_ConsensusSeq(filename=f'{samplename}_consensus_sequence.fasta', genome_consensus=genome_consensus['unmasked'])
		write_ConsensusSeq(filename=f'{samplename}_masked_consensus_sequence.fasta', genome_consensus=genome_consensus['masked'])

	remove_IntermediateFiles()

	logger.logger.info('Finished running detect_variants.py\n')







