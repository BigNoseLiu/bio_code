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
	with open(infile) as inf, open(outfile, 'w') as outf:
		outf.write('tax_id\trefname\tgene\tVFgene\tsummary\tspecies\tidentity\ttype_judge\n')
		inf.readline()
		for l in inf.readlines():
			l = l.strip().split('\t')
			if l[0] in Anno and float(l[6]) >= identity:
				Info[l[0]] = float(l[6])
	
		for i in sorted(Info.items(), key=lambda x: (-x[1], x[0])):
			outf.write(f'{Anno[i[0]]}\t{Info[i[0]]}\t-\n')

if __name__ == '__main__':
	args = readArg(argv)
	infile = os.path.abspath(args['in'])
	outfile = os.path.abspath(args['out'])

	filterResult(infile, outfile, args['identity'])
