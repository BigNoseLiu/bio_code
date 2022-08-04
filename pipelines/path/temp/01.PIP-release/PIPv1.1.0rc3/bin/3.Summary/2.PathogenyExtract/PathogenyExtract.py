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
getSpeciesReads = os.path.join(Bin, 'getSpeciesReads.pl')
#####################################################

def readArg(args):
	parser = ArgumentParser(description = \
		f'{os.path.basename(os.path.abspath(__file__))} {__version__} ({__date__}) {__author__} {__email__}')
	args = parser.add_argument
	args('--list', type=str, action='store', required=True, \
		help='Sample information table, format: "sample\\tDNA|RNA\\tread1[\\tread2]\\tread_length\\tkout \
        \\tkreport\\tmpa.kreport\\ttaxid.list\\tTotal_Detail.txt\\tPathogeny_Detail.txt".')
	args('--index', type=str, action='store', default='RPM', choices=['RPM', 'NV'], \
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
	with open(L) as inf, open(SH, 'w') as sh:
		for l in inf.readlines():
			line += 1
			l = l.strip().split('\t')
			samplePath = outdir / l[0] / l[1]
			createPath(samplePath)

			if Type == 'SE':
				sh.write(f'{getSpeciesReads} --fq {l[2]} --taxid {l[-3]} --krakenout {l[-6]} --outf fastq --outdir {samplePath}\n')
			else:
				sh.write(f'{getSpeciesReads} --fq {l[2]},{l[3]} --taxid {l[-3]} --krakenout {l[-6]} --outf fastq --outdir {samplePath}\n')

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
		os.system(f'{cfg.sbatch} --maxjob 30 --jobprefix PathogenyExtract {cmd}')
	else:
		os.system(f'{cfg.watchDog} --mem 1G --num_paral 30 {cmd}')
	printInfo(os.path.abspath(__file__), 2)
