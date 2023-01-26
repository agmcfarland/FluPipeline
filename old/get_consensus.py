def write_ConsensusSeq(filename, genome_consensus):
	'''
	'''
	with open(filename, 'w') as outfile:
		for chrom, seq in genome_consensus.items():
			outfile.write(f'>{chrom}\n')
			[outfile.write(i) if e%60!=0 else outfile.write(f'{i}\n') for e, i in enumerate(seq, 1)]
			outfile.write('\n')

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