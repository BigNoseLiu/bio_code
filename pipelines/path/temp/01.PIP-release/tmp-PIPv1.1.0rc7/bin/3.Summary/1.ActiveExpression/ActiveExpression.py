#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'qiuwancen'
__email__ = '972538446@qq.com'
__date__ = '2021-06-16'
__version__ = 'V1.0'

import os
import sys
from sys import argv
from pathlib import Path
from argparse import ArgumentParser

Bin = os.path.dirname(os.path.abspath(__file__))

sys.path.append(os.path.join(Bin, '../../../conf'))
import configure as cfg
from configure import createPath, checkFile, checkType, printInfo

###################### scripts ######################
ActiveExpressionStat = os.path.join(Bin, 'ActiveExpressionStat.pl')
#####################################################

def readArg(args):
	parser = ArgumentParser(description = \
		f'{os.path.basename(os.path.abspath(__file__))} {__version__} ({__date__}) {__email__}')
	args = parser.add_argument
	args('--list', type=str, action='store', required=True, \
		help='Sample information table, format: "sample\\tDNA|RNA\\tread1[\\tread2]\\tread_length\\tkout \
        \\tkreport\\tmpa.kreport\\ttaxid.list\\tTotal_Detail.txt\\tPathogeny_Detail.txt".')
	args('--index', type=str, action='store', default='RPM', choices=['RPM', 'TPM'], \
		help='Index used for uniformization. [%(default)s]')
	args('--thread', type=int, action='store', default=1, \
		help='Thread. [%(default)s]')
	args('--type', type=str, action='store', default='SE', choices=['SE', 'PE'], \
		help='Input sequence type. [%(default)s]')
	args('--run_mode', type=str, action='store', default='local', choices=['local', 'Slurm'], \
		help='Set the run mode, local or Slurm, the option --partition must be set if you set run mode to Slurm. [%(default)s]')
	args('--partition', type=str, action='store', choices=['ngs-node1', 'ngs-node2', 'ngs-node3'], \
		help='partition requested, it need to be set if you set --run_mode to Slurm.')
	args('--outdir', type=str, action='store', default='.', \
		help='Output path. [%(default)s]')
	args('--version', action='version', version=f'{os.path.basename(os.path.abspath(__file__))} {__version__} ({__date__})', \
		help='Print the current version and exit.')
	return vars(parser.parse_args())

def generateScripts(L, outdir):
	outdir = Path(outdir)
	shell = outdir / 'shell'
	createPath(outdir, shell)

	SH = shell / 'ActiveExpression.sh'
	line = 0
	Link = {}
	with open(L) as inf:
		for l in inf.readlines():
			l = l.strip().split('\t')

			if l[0] in Link:
				Link[l[0]].update({l[1]: {'Total': l[-2], 'Pathogeny': l[-1]}})
			else:
				Link.update({l[0]: {l[1]: {'Total': l[-2], 'Pathogeny': l[-1]}}})

	with open(SH, 'w') as sh:
		for sample in Link:
			if 'DNA' in Link[sample] and 'RNA' in Link[sample]:
				line += 1
				samplePath = outdir / sample
				createPath(samplePath)

				sh.write('%(perl)s %(ActiveExpressionStat)s --DNA %(TotalDNA)s --RNA %(TotalRNA)s '
				'--index %(index)s --out %(samplePath)s/Total.ActiveExpression.txt\n'
				'%(perl)s %(ActiveExpressionStat)s --DNA %(PathogenyDNA)s --RNA %(PathogenyRNA)s '
				'--index %(index)s --out %(samplePath)s/Pathogeny.ActiveExpression.txt\n' % \
				{'perl': cfg.perl, 'ActiveExpressionStat': ActiveExpressionStat, 'TotalDNA': Link[sample]["DNA"]["Total"],
				'TotalRNA': Link[sample]["RNA"]["Total"], 'PathogenyDNA': Link[sample]["DNA"]["Pathogeny"], 
				'PathogenyRNA': Link[sample]["RNA"]["Pathogeny"], 'index': index, 'samplePath': samplePath})

	return SH, line

if __name__ == '__main__':
	args = readArg(argv)
	infile = os.path.abspath(args['list'])
	outdir = os.path.abspath(args['outdir'])
	Type = args['type']
	mode = args['run_mode']
	index = args['index']

	checkFile(infile)

#	if mode.upper() == 'SLURM' and not args['partition']:
#		print('Error: the option --partition must be set when you set --run_mode to Slurm')
#		exit(0)

	cmd, lineNum = generateScripts(infile, outdir)
	printInfo(os.path.abspath(__file__), 1)
	if mode.upper() == 'SLURM':
		os.system(f'{cfg.perl} {cfg.sbatch} --maxjob 140 --jobprefix ActiveExpression {cmd}')
	else:
		os.system(f'{cfg.perl} {cfg.watchDog} --mem 1G --num_paral 140 {cmd}')
	printInfo(os.path.abspath(__file__), 2)
