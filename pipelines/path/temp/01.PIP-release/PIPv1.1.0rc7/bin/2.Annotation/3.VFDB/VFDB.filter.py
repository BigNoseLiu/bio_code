#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'qiuwancen'
__email__ = '972538446@qq.com'
__date__ = '2021-07-20'
__version__ = 'V1.0'

import os
import sys
from sys import argv
from argparse import ArgumentParser
import pandas as pd
Bin = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(Bin, '../../../conf'))
import configure as cfg

def readArg(args):
	parser = ArgumentParser(description = \
		f'{os.path.basename(os.path.abspath(__file__))} {__version__} ({__date__}) {__email__}')
	args = parser.add_argument
	args('--in', type=str, action='store', required=True, \
		help='Enter the results of virulence gene analysis.')
	args('--db', type=str, action='store', default='20210709', \
		help='Database version of VFDB. [%(default)s]')
	args('--identity', type=float, action='store', default=95.0, \
		help='Identity threshold. [%(default)s]')
	args('--out', type=str, action='store', default='./out.txt', \
		help='Output file. [%(default)s]')
	args('--version', action='version', version=f'{os.path.basename(os.path.abspath(__file__))} {__version__} ({__date__})', \
		help='Print the current version and exit.')
	return vars(parser.parse_args())

def filterResult(infile, outfile, identity=95.0):
	Anno = {}
	with open(f'{cfg.VFDB}/{args["db"]}/VFDB.Anno.txt') as anno:
		anno.readline()
		for line in anno.readlines():
			line = line.strip()
			lines = line.split('\t')
			Anno[lines[1]] = line
	Info = {}
	mapstat={}
	with open(infile) as inf,open(infile.split('.res')[0]+'.mapstat') as inf2,open(outfile, 'w') as outf:
		outf.write('Gene\tReads_Num\tCoverage(%)\tIdentity(%)\tPathogen\tDatabaseID\tAnnotation\n')
		inf.readline()
		inf2.readline()
		for l2 in inf2.readlines()[6:]:
			l2 = l2.strip().split('\t')
			mapstat[l2[0]] = float(l2[1])
		for l in inf.readlines():
			l = l.strip().split('\t')
			if l[0] in Anno and float(l[6]) >= identity:
				Info[l[0]] = float(l[6])
				gene=Anno[l[0]].split('\t')[2]
				species=Anno[l[0]].split('\t')[-1]
				annt=Anno[l[0]].split('\t')[-2]
				outf.write(f"{gene}\t{mapstat[l[0]]}\t{l[5]}\t{l[6]}\t{species}\t{l[0]}\t{annt}\n")

if __name__ == '__main__':
	args = readArg(argv)
	infile = os.path.abspath(args['in'])
	outfile = os.path.abspath(args['out'])
	filterResult(infile, outfile, args['identity'])
	from pathlib import Path
	pd.read_table(outfile).to_excel(Path(outfile).with_suffix('.xlsx'), index=False)
