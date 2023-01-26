import sys
import os 
from os.path import join as pjoin
import multiprocessing as mp
from multiprocessing import Pool
import pandas as pd
import re
import glob
import shutil
import subprocess
from datetime import datetime
from Bio.Seq import Seq
from Bio import SeqIO, SeqFeature
import argparse
import logging
import numpy as np
from .processing_classes import SequencingSample, RunLogger


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


def reference_NextCladeLookUp(reference):
	'''
	Returns the the nextclade folder name containing all necessary comparison data to run nextclade.
	If the reference is not available in the dictionary will return none available. this is useful when user is 
	supplying their own references.
	'''
	reference_dict = {
	'IBV_Victoria_ref':'flu_vic_ha',
	'IBV_Yamagata_ref':'flu_yam_ha',
	'H1N1pdm_ref':'flu_h1n1pdm_ha',
	'H3N2_ref':'flu_h3n2_ha',
	'H6N2_Massachusetts':'none available',
	'H1N1_A_Brisbane_59_2007':'none available'
	}

	try:
		refname = reference_dict[reference]
	except:
		return('none available')

	return(refname)


def select_BestReference(samplename, directory=os.getcwd()):
	'''
	Selects the best reference from the output of summarize_ReferenceCoverage
	(most segments, highest coverage, highest read depth). 
	'''
	df = pd.read_csv(pjoin(directory,'reference_ha_coverage_{}.csv'.format(samplename)))
	df = df.sort_values(['segment_coverage','average_read_depth'], ascending = [False,False])
	return(df.head(1)['reference'].item())


if __name__=='__main__':
	pass
















