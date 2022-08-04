#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'qiuwancen'
__email__ = '972538446@qq.com'
__date__ = '2021-12-04'
__version__ = 'V1.0'

import os
import sys
import math
import pandas as pd
from sys import argv
from pathlib import Path
Bin = Path(__file__).resolve().parent
sys.path.append(os.path.join(Bin, '../../conf'))
import configure as cfg
from configure import createPath
from argparse import ArgumentParser

def readArg(args):
	parser = ArgumentParser(description='{} {} ({}) {}'.format(\
		os.path.basename(os.path.abspath(__file__)), __version__, __date__, __email__))
	args = parser.add_argument
	args('--qc', type=str, action='store', required=True, \
		help='QC table.')
	args('--anno', type=str, action='store', required=True, \
		help='List of annotation.')
	args('--volume', type=int, action='store', default=10, \
		help='The reads numer that normalize to (Unit: M). [%(default)s]')
	args('--outdir', type=str, action='store', default='./', \
		help='Output Path. [%(default)s]')
	args('--version', action='version', version='{} {} ({})'.format(\
		os.path.basename(os.path.abspath(__file__)), __version__, __date__), \
		help='Print the current version and exit.')
	return vars(parser.parse_args())

def dealTable(infile, fold):
	out_cols = ['Kingdom', 'GenusName', 'GenusReads', 'GroupName', 'GroupReads', 'ScientificName', 	'ChineseName',
				'Reads_Number', 'GMRN', 'TPM', 'RH-RPM','Normalization', 'Verification', 'Focus']
	
	table = pd.read_table(infile, header=0, delimiter='\t')
	# code block for replacing the next code line
	if not table.empty:
		table['Normalization'] = table.apply(lambda x: math.ceil(int(x['Reads_Number']) * fold), axis=1)
	else:
		table = pd.DataFrame(columns=out_cols)
	table.to_csv(f'{infile}.2',sep='\t', index=False, columns=out_cols)		
	# table.insert(11, 'Normalization', [math.ceil(int(i)*fold) for i in list(table['Reads_Number'])])
	
	os.system(f'mv {infile}.2 {infile}')

if __name__ == '__main__':
	args = readArg(argv)
	qc = os.path.abspath(args['qc'])
	anno = os.path.abspath(args['anno'])
	outdir = os.path.abspath(args['outdir'])
	volume = args['volume'] * 1000000
	createPath(outdir)

	# Statistical sample data volume
	df = pd.read_table(qc, header=0, index_col='Sample', delimiter="\t")
	Rows = list(df.index)
	samples = list(df.columns)
	types = list(df.iloc[Rows.index('Type')])
	sampleType = ['_'.join(i) for i in list(zip(samples, types))]
	RawReads = list(df.iloc[Rows.index('Raw Reads')])
	Sample2Reads = dict(zip(sampleType, RawReads))

	# Normalize the sample reads number to N (M)
	with open(anno) as List:
		for line in List.readlines():
			line = line.strip().split('\t')
			times = volume / int(Sample2Reads[f'{line[0]}_{line[1]}'])
			dealTable(line[-1], times)
			dealTable(line[-2], times)


