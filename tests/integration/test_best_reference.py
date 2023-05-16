"""
Test that the best_reference module selects the best Flu reference genome.
"""

import unittest
import tempfile
import os
from os.path import join as pjoin
import subprocess
import shutil
import pandas as pd
import sys
sys.path.append('../../')
from pylibs.processing_functions import select_BestReference

class TestBestReference(unittest.TestCase):
	"""
	Test that the best_reference module selects the best Flu reference genome.
	"""
	def setUp(self):
		"""
		Set up common paths for tests.
		"""
		self.flu_pipeline_directory = os.path.split(os.path.realpath(os.path.abspath(__file__)))[0].replace('/tests/integration','')
		self.tests_data = pjoin(self.flu_pipeline_directory, 'tests', 'data')
		self.reference_strains_dir = pjoin(self.flu_pipeline_directory, 'references')

		print(f'FluPipeline dirctory: {self.flu_pipeline_directory}')

	# class self():
		"""
		Interactive troubleshooting
		"""
	# 	def __init__(self):
	# 		self.flu_pipeline_directory = '/data/FluPipeline'
	# 		self.tests_data = pjoin(self.flu_pipeline_directory, 'tests', 'data')
	# 		self.reference_strains_dir = pjoin(self.flu_pipeline_directory, 'references')
	# self = self()

	def test_best_reference_module(self):
		"""
		Test that the module best_reference has the correct number of rows and picks the best reference for the test sample.
		"""
		self.output_dir = tempfile.mkdtemp(dir = tempfile.gettempdir())
		sample_name = os.path.basename(self.output_dir)
		r1_path = pjoin(self.tests_data, 'ds_filt_test1_R1_001.fastq.gz')
		r2_path = pjoin(self.tests_data, 'ds_filt_test1_R2_001.fastq.gz')

		subprocess.run(['python', '-m', 'bin.best_reference',
			'--baseDir', self.output_dir,
			'--logDir', self.output_dir,
			'--referenceStrainsDir', self.reference_strains_dir,
			'--R1', r1_path,
			'--R2', r2_path])

		self.assertEqual(
			pd.read_csv(pjoin(self.output_dir, f'combined_coverage_stats_{sample_name}.csv')).shape, (38, 4))

		self.assertEqual(
			pd.read_csv(pjoin(self.output_dir, f'reference_coverage_{sample_name}.csv')).shape, (5, 5))

		self.assertEqual(
			pd.read_csv(pjoin(self.output_dir, f'reference_ha_coverage_{sample_name}.csv')).shape, (5, 5))

		self.assertEqual(
			select_BestReference(samplename = sample_name, directory = self.output_dir), 'H1N1pdm_ref.fasta')


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

