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
from configure import checkFile, printInfo, checkType, createPath

def readArg(args):
	parser = ArgumentParser(description = \
		f'{os.path.basename(os.path.abspath(__file__))} {__version__} ({__date__}) {__email__}')
	args = parser.add_argument
	args('--list', type=str, action='store', required=True, \
		help='List of Fastq, format: "sample\\tDNA|RNA\\tread1[\\tread2]\\tread_length".')
	args('--thread', type=int, action='store', default=8, \
		help='Thread. [%(default)s]')
	args('--threshold', type=float, action='store', default=0.3, \
		help='Coverage threshold. [%(default)s]')
	args('--identity', type=float, action='store', default=0.6, \
		help='Identity. [%(default)s]')
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

def generateScript(L, outdir):
	outdir = Path(outdir)
	shell = outdir / 'shell'
	createPath(outdir, shell)

	SH = shell / 'ResistanceGene.sh'
	sList = outdir / 'sample.list'
	line = 0

	envs = os.path.dirname(cfg.freebayesP)
	with open(L) as inf, open(sList, 'w') as sl, open(SH, 'w') as sh:
		for l in inf.readlines():
			line += 1
			l = l.strip().split('\t')
			checkType(Type, l)

			samplePath = outdir / l[0] / l[1]
			createPath(samplePath)

			if Type == 'SE':
				sh.write('%(python39)s %(resfinder)s --cpu %(thread)s -db_res_kma %(resfinderDB)s -s Other '
					'-l %(coverage)s -t %(identity)s --acquired -o %(samplePath)s -ifq %(fastq)s\n' % \
					{'resfinder': cfg.resfinder, 'resfinderDB': cfg.resfinderDB, 'thread': args['thread'], 'fastq': l[2], \
					'coverage': args['threshold'], 'identity': args['identity'], 'samplePath': samplePath, 'python39': cfg.python39})
			else:
				sh.write('%(python39)s %(resfinder)s --cpu %(thread)s -db_res_kma %(resfinderDB)s -s Other '
					'-l %(coverage)s -t %(identity)s --acquired -o %(samplePath)s -ifq %(fastq1)s %(fastq2)s\n' % \
					{'resfinder': cfg.resfinder, 'resfinderDB': cfg.resfinderDB, 'thread': args['thread'], 'fastq1': l[2], 'fastq2': l[3], \
					'coverage': args['threshold'], 'identity': args['identity'], 'samplePath': samplePath, 'python39': cfg.python39})
			sh.write('mv %(samplePath)s/ResFinder_results_tab.txt %(samplePath)s/ARG.cov.pick.xls &\n'
				'rm -fr %(samplePath)s/resfinder_kma &\n'
				'wait\n' % \
				{'samplePath': samplePath})
	return SH, line

if __name__ == '__main__':
	args = readArg(argv)
	infile = os.path.abspath(args['list'])
	outdir = os.path.abspath(args['outdir'])
	Type = args['type']
	mode = args['run_mode']

	checkFile(infile)

#	if mode.upper() == 'SLURM' and not args['partition']:
#		print('Error: the option --partition must be set when you set --run_mode to Slurm')
#		exit(0)

	cmd, lineNum = generateScript(infile, outdir)
	printInfo(os.path.abspath(__file__), 1)
	if mode.upper() == 'SLURM':
		os.system(f'{cfg.perl} {cfg.sbatch} --maxjob 51 --cpu {args["thread"]} --lines 4 --jobprefix ResistanceGene {cmd}')
	else:
		os.system(f'{cfg.perl} {cfg.watchDog} --mem 2G --num_paral 17 --num_line 4 {cmd}')
	printInfo(os.path.abspath(__file__), 2)
