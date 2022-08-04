#!/usr/bin/env python
# -* coding: utf-8 -*-

__author__ = 'qiuwancen'
__email__ = '972538446@qq.com'
__date__ = '2021-01-01'
__version__ = 'V1.0'

import os
import sys
import fileinput
from sys import argv
from argparse import ArgumentParser
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../conf'))
import configure as cfg
from configure import checkFile

def readArg(args):
	parser = ArgumentParser(description = \
		f'{os.path.basename(os.path.abspath(__file__))} {__version__} ({__date__}) {__author__} {__email__}')
	args = parser.add_argument
	args('--ko', type=str, action='store', required=True, \
		help="Input kraken's output.")
	args('--db', type=str, action='store', default='PIDB202101', \
		help='Database version of PIP. [%(default)s]')
	args('--fc', type=str, action='store', default='min', choices=['min', 'highconf', 'unique'], \
		help='Filter criteria. [%(default)s]')
	args('--out', type=str, action='store', default='./out.txt', \
		help='Output. [%(default)s]')
	args('--version', action='version', version=f'{os.path.basename(os.path.abspath(__file__))} {__version__} ({__date__})', \
        help='Print the current version and exit.')
	return vars(parser.parse_args())

def filterOut(ko, cf, out, criteria='min'):
#	Threshold = {}
	Species = {}
	with open(cf + '/seqTaxid2speciesTaxid') as idInfo:
		for idl in idInfo.readlines():
			idl = idl.strip().split('\t')
			if idl[1] != 'NULL':
				Species[idl[0]] = idl[1]

	with open(cf + '/all.Taxonomy.txt') as cutoff:
		cutoff.readline()
		for cl in cutoff.readlines():
			cl = cl.strip().split('\t')
			Species[cl[0]] = cl[1]
#			if cl[0] not in Threshold:
#				Threshold.update({cl[0]: {'min': 0, 'highconf': 0}})
#			Threshold[cl[0]]['min'] = int(cl[2])
#			Threshold[cl[0]]['highconf'] = int(cl[3])

	with open(out, 'w') as outf:
		for kl in fileinput.input([ko]):
			if kl.startswith('U'):
				outf.write(kl)
			elif criteria == 'unique':
				SID = []
				kL = kl.strip().split('\t')
				if kL[2] not in Species:
					outf.write(kl)
				else:
					for detail in kL[-1].split(' '):
						key, value = detail.split(':')
						if key != '|' and key in Species:
							SID.append(Species[key])
					SID = list(set(SID))
					if len(SID) > 1:
						outf.write(f'C\t{kL[1]}\t131567\t{kL[3]}\t{kL[4]}\n')
					else:
						outf.write(kl)
#			else:
#				num = flag = 0
#				kL = kl.strip().split('\t')
#				for detail in kL[-1].split(' '):
#					key, value = detail.split(':')
#					if key != '|':
#						if key == kL[2]:
#							num += int(value)
#							flag = 1
#				if flag == 0:
#					outf.write(kl)
#					continue
#				if kL[2] not in Threshold:
#					outf.write(kl)
#					continue
#				if kL[2] in Threshold:
#					if num >= Threshold[kL[2]][criteria]:
#						outf.write(kl)
#					else:
#						outf.write(f'C\t{kL[1]}\t131567\t{kL[3]}\t{kL[4]}\n')

if __name__ == '__main__':
	args = readArg(argv)
	kout = os.path.abspath(args['ko'])
	outfile = os.path.abspath(args['out'])

	checkFile(kout)
	
	if args['fc'] not in ['min', 'highconf', 'unique']:
		print("There are problems with filtering criteria.")
		exit(1)
	dbFilePath = os.path.join(cfg.KrakenDB, args['db'])
#	cutFile = os.path.abspath(args['db'])
	filterOut(kout, dbFilePath, outfile, args['fc'])
