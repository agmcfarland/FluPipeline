"""
Test that the detect_variants module calls the expected number of variants.
"""

import unittest
import tempfile
import os
from os.path import join as pjoin
import subprocess
import shutil
import pandas as pd
import sys
from Bio import SeqIO
from Bio.Seq import Seq
sys.path.append('../../')
from pylibs.processing_functions import convert_GBKToFasta

class TestDetectVariants(unittest.TestCase):
	"""
	Check that the detect_variants modules calls the expected number of variants.
	"""

	def setUp(self):
		"""
		Set up common paths for tests.
		"""
		self.flu_pipeline_directory = os.path.split(os.path.realpath(os.path.abspath(__file__)))[0].replace('/tests/integration','')
		self.tests_data = pjoin(self.flu_pipeline_directory, 'tests', 'data')
		self.R1 = pjoin(self.tests_data, 'ds_filt_test1_R1_001.fastq.gz')
		self.R2 = pjoin(self.tests_data, 'ds_filt_test1_R2_001.fastq.gz')
		self.reference_genome = pjoin(self.flu_pipeline_directory, 'references', 'H1N1pdm_ref')
		convert_GBKToFasta(filename = self.reference_genome)

		with open(self.reference_genome + '.fasta') as input_handle:
			for record in SeqIO.parse(input_handle, 'fasta'):
				if record.id == 'CY266198_1_1':
					self.reference_record = record
					break

	def test_convert_gbk_to_fasta(self):
		"""
		Test that convert_gbk_to_fasta converts a GBK file to a fasta file with sequences of expected names and sizes.
		"""
		expected_data = {
			'CY266198_1_1' : 2314,
			'CY266197_1_2' : 2296,
			'CY266196_1_3' : 2175,
			'CY266195_1_8' : 860,
			'CY266194_1_5' : 1541,
			'CY266193_1_6' : 1420,
			'CY266192_1_7' : 992,
			'CY266191_1_4' : 1734
		}
		for record in SeqIO.parse(self.reference_genome + '.fasta', 'fasta'):
			seq_length = len(str(record.seq))
			self.assertEqual(
				expected_data[record.id], seq_length)

	def test_detect_variants_bcftools(self):
		"""
		Test bcftools variant caller. Check that the variants called at different filter thresholds are
		consistent. Check that for segment 1 the bases in the consensus sequence match the majorVariants
		table.
		"""
		print(os.path.exists(self.R1))
		self.output_dir = tempfile.mkdtemp(dir = tempfile.gettempdir())
		sample_name = os.path.basename(self.output_dir)
		print(self.output_dir)

		subprocess.run((
			'python', '-m', 'bin.detect_variants',
			'--baseDir', self.output_dir,
			'--logDir', self.output_dir,
			'--variant_caller', 'bcftools',
			'--refGenomeFasta', self.reference_genome + '.fasta',
			'--R1', self.R1,
			'--R2', self.R2,
			'--consensus_sequence'))


		self.assertEqual(pd.read_csv(
			pjoin(self.output_dir, sample_name + '_allVariants.csv')).shape, (324,30))

		self.assertEqual(pd.read_csv(
			pjoin(self.output_dir, sample_name + '_passingVariants.csv')).shape, (310,30))

		self.assertEqual(pd.read_csv(
			pjoin(self.output_dir, sample_name + '_majorVariants.csv')).shape, (310,30))

		for record in SeqIO.parse(self.reference_genome + '.fasta', 'fasta'):
			if record.id == 'CY266198_1_1':
				reference_record = record
				break

		with open(pjoin(self.output_dir, sample_name + '_consensus_sequence.fasta')) as input_handle:
			for record in SeqIO.parse(input_handle, 'fasta'):
				if record.id == 'CY266198_1_1':
					variant_record = record
					break

		mismatch_dict = {}
		str_variant_record = str(variant_record.seq)
		for e, base in enumerate(str(reference_record.seq)):
			variant_base = str_variant_record[e]
			if variant_base != base:
				mismatch_dict[e + 1] = [base, variant_base]

		df_majorVariants = pd.read_csv(pjoin(self.output_dir, sample_name + '_majorVariants.csv'))
		print('Major variants called:', len(df_majorVariants))
		
		df_majorVariants = df_majorVariants[df_majorVariants['CHROM'] == 'CY266198_1_1']
		
		for e, row in df_majorVariants.iterrows():
			mismatch_data = mismatch_dict[row['POS']]
			self.assertEqual(row['REF'], mismatch_data[0])
			self.assertEqual(row['ALT'], mismatch_data[1])

	def test_detect_variants_bbtools(self):
		"""
		Test bbtools variant caller. Check that the variants called at different filter thresholds are
		consistent. Check that for segment 1 the bases in the consensus sequence match the majorVariants
		table.
		"""
		print(os.path.exists(self.R1))
		self.output_dir = tempfile.mkdtemp(dir = tempfile.gettempdir())
		sample_name = os.path.basename(self.output_dir)
		print(self.output_dir)

		subprocess.run((
			'python', '-m', 'bin.detect_variants',
			'--baseDir', self.output_dir,
			'--logDir', self.output_dir,
			'--variant_caller', 'bbtools',
			'--refGenomeFasta', self.reference_genome + '.fasta',
			'--R1', self.R1,
			'--R2', self.R2,
			'--consensus_sequence'))

		self.assertEqual(pd.read_csv(
			pjoin(self.output_dir, sample_name + '_allVariants.csv')).shape, (1757,37))

		self.assertEqual(pd.read_csv(
			pjoin(self.output_dir, sample_name + '_passingVariants.csv')).shape, (316,37))

		self.assertEqual(pd.read_csv(
			pjoin(self.output_dir, sample_name + '_majorVariants.csv')).shape, (308,37))


		with open(pjoin(self.output_dir, sample_name + '_consensus_sequence.fasta')) as input_handle:
			for record in SeqIO.parse(input_handle, 'fasta'):
				if record.id == 'CY266198_1_1':
					variant_record = record
					break

		mismatch_dict = {}
		str_variant_record = str(variant_record.seq)
		for e, base in enumerate(str(self.reference_record.seq)):
			variant_base = str_variant_record[e]
			if variant_base != base:
				mismatch_dict[e + 1] = [base, variant_base]

		df_majorVariants = pd.read_csv(pjoin(self.output_dir, sample_name + '_majorVariants.csv'))
		print('Major variants called:', len(df_majorVariants))
		
		df_majorVariants = df_majorVariants[df_majorVariants['CHROM'] == 'CY266198_1_1']
		
		for e, row in df_majorVariants.iterrows():
			mismatch_data = mismatch_dict[row['POS']]
			self.assertEqual(row['REF'], mismatch_data[0])
			self.assertEqual(row['ALT'], mismatch_data[1])

	def test_detect_variants_lofreq(self):
		"""
		Test lofreq variant caller. Check that the variants called at different filter thresholds are
		consistent. Check that for segment 1 the bases in the consensus sequence match the majorVariants
		table.
		"""
		print(os.path.exists(self.R1))
		self.output_dir = tempfile.mkdtemp(dir = tempfile.gettempdir())
		sample_name = os.path.basename(self.output_dir)
		print(self.output_dir)

		subprocess.run((
			'python', '-m', 'bin.detect_variants',
			'--baseDir', self.output_dir,
			'--logDir', self.output_dir,
			'--variant_caller', 'lofreq',
			'--refGenomeFasta', self.reference_genome + '.fasta',
			'--R1', self.R1,
			'--R2', self.R2,
			'--consensus_sequence'))

		self.assertEqual(pd.read_csv(
			pjoin(self.output_dir, sample_name + '_allVariants.csv')).shape, (317,14))

		self.assertEqual(pd.read_csv(
			pjoin(self.output_dir, sample_name + '_passingVariants.csv')).shape, (309,14))

		self.assertEqual(pd.read_csv(
			pjoin(self.output_dir, sample_name + '_majorVariants.csv')).shape, (308,14))

		with open(pjoin(self.output_dir, sample_name + '_consensus_sequence.fasta')) as input_handle:
			for record in SeqIO.parse(input_handle, 'fasta'):
				if record.id == 'CY266198_1_1':
					variant_record = record
					break

		mismatch_dict = {}
		str_variant_record = str(variant_record.seq)
		for e, base in enumerate(str(self.reference_record.seq)):
			variant_base = str_variant_record[e]
			if variant_base != base:
				mismatch_dict[e + 1] = [base, variant_base]

		df_majorVariants = pd.read_csv(pjoin(self.output_dir, sample_name + '_majorVariants.csv'))
		print('Major variants called:', len(df_majorVariants))
		
		df_majorVariants = df_majorVariants[df_majorVariants['CHROM'] == 'CY266198_1_1']
		
		for e, row in df_majorVariants.iterrows():
			mismatch_data = mismatch_dict[row['POS']]
			self.assertEqual(row['REF'], mismatch_data[0])
			self.assertEqual(row['ALT'], mismatch_data[1])

	def test_detect_variants_freebayes(self):
		"""
		Test freebayes variant caller. Check that the variants called at different filter thresholds are
		consistent. Check that for segment 1 the bases in the consensus sequence match the majorVariants
		table.

		This test is incomplete.
		"""
		print(os.path.exists(self.R1))
		self.output_dir = tempfile.mkdtemp(dir = tempfile.gettempdir())
		sample_name = os.path.basename(self.output_dir)
		print(self.output_dir)

		subprocess.run((
			'python', '-m', 'bin.detect_variants',
			'--baseDir', self.output_dir,
			'--logDir', self.output_dir,
			'--variant_caller', 'freebayes',
			'--refGenomeFasta', self.reference_genome + '.fasta',
			'--R1', self.R1,
			'--R2', self.R2,
			'--consensus_sequence'))


		self.assertEqual(pd.read_csv(
			pjoin(self.output_dir, sample_name + '_allVariants.csv')).shape, (301,54))

		self.assertEqual(pd.read_csv(
			pjoin(self.output_dir, sample_name + '_passingVariants.csv')).shape, (288,54))

		self.assertEqual(pd.read_csv(
			pjoin(self.output_dir, sample_name + '_majorVariants.csv')).shape, (285,54))


	def tearDown(self):
		"""
		Remove test data output files and directories that were created.
		"""
		try:
			shutil.rmtree(self.output_dir)
		except:
			print('Could not delete temporary directory')



if __name__ == '__main__':
	unittest.main()
