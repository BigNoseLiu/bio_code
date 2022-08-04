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

Bin = Path(__file__).resolve().parent
sys.path.append(os.path.join(Bin, '../../conf'))
import configure as cfg
from configure import checkFile, createPath, printInfo

def readArg(args):
	parser = ArgumentParser(description = f'{Path(__file__).name} {__version__} ({__date__}) {__email__}')
	args = parser.add_argument
	args('--list', type=str, action='store', required=True, \
		help='List of Fastq, format: "sample\\tDNA|RNA\\tread1[\\tread2]\\tread_length".')
	args('--analysis', type=str, action='store', default='Fastp,RemoveHost,rRNAFilter', \
		help='QC types including Fastp, RemoveHost that separated by ",". [%(default)s]')
	args('--type', type=str, action='store', default='SE', choices=['SE', 'PE'], \
		help='Input sequence type. [%(default)s]')
	args('--run_mode', type=str, action='store', default='local', choices=['local', 'Slurm', 'SGE'], \
		help='Set the run mode, local, Slurm or SGE, the option --partition must be set if you set run mode to Slurm or SGE. [%(default)s]')
	args('--partition', type=str, action='store', choices=['ngs-node1', 'ngs-node2', 'ngs-node3'], \
		help='partition requested, it need to be set if you set --run_mode to Slurm or SGE.')
	args('--outdir', type=str, action='store', default='.', \
		help='Output path. [%(default)s]')
	args('--version', action='version', version=f'{Path(__file__).name} {__version__} ({__date__})', \
		help='Print the current version and exit.')
	return vars(parser.parse_args())

def generateScripts(L, outdir, analysis):
	Scripts = {
		'Fastp': Bin / '1.Fastp/Fastp.py', 
		'RemoveHost': Bin / '2.RemoveHost/RemoveHost.py',
		'rRNAFilter': Bin / '3.rRNAFilter/rRNAFilter.py',
	}
	Workdir = {
		'Fastp': outdir / '01.Fastp',
		'RemoveHost': outdir / '02.RemoveHost',
		'rRNAFilter': outdir / '03.rRNAFilter',
	}
	shell = outdir / 'shell'
	createPath(outdir, shell)
	with open(shell / 'QC.sh', 'w') as sh:
		for step in analysis.split(','):
			sh.write(f'{cfg.python} {Scripts[step]} --list {L} --type {Type} {mode} --outdir {Workdir[step]}\n')
			L = str(Workdir[step]) + '/sample.list'
	return str(shell) + '/QC.sh', L

if __name__ == '__main__':
	args = readArg(argv)
	infile = Path(args['list']).resolve()
	analysis = args['analysis']
	Type = args['type']
	mode = args['run_mode']
	outdir = Path(args['outdir']).resolve()

	checkFile(infile)

#	if mode.upper() == 'SLURM' and not args['partition']:
#		print('Error: the option --partition must be set when you set --run_mode to Slurm')
#		exit(0)
	if mode.upper() == 'SLURM' or mode.upper() == 'SGE':
		mode = f"--run_mode {mode}"
		if args['partition']:
			mode += f" --partition {args['partition']}"
	else:
		mode = f"--run_mode {mode}"

	cmd, list = generateScripts(infile, outdir, analysis)
#	printInfo(os.path.abspath(__file__), 1)
	os.system(f'sh {cmd}')
	os.system(f'ln -s -f {list} {outdir}/sample.list')
#	printInfo(os.path.abspath(__file__), 2)
