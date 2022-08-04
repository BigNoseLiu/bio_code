#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'qiuwancen'
__email__ = '972538446@qq.com'
__date__ = '2021-12-01'
__version__ = 'V1.0'

import os
import re
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
	args('--analysis', type=str, action='store', default='Taxonomy=PIDB202101_Middle,ResistanceGene,VFDB=20210709', \
		help='Annotation types including Taxonomy, ResistanceGene that separated by ",", the \
		database version of Taxonomy follow the "=". [%(default)s]')
	args('--type', type=str, action='store', default='SE', choices=['SE', 'PE'], \
		help='Input sequence type. [%(default)s]')
	args('--run_mode', type=str, action='store', default='local', choices=['local', 'Slurm'], \
		help='Set the run mode, local or Slurm, the option --partition must be set if you set run mode to Slurm. [%(default)s]')
	args('--partition', type=str, action='store', choices=['ngs-node1', 'ngs-node2', 'ngs-node3'], \
		help='partition requested, it need to be set if you set --run_mode to Slurm.')
	args('--outdir', type=str, action='store', default='.', \
		help='Output path. [%(default)s]')
	args('--version', action='version', version=f'{Path(__file__).name} {__version__} ({__date__})', \
		help='Print the current version and exit.')
	return vars(parser.parse_args())

def generateScripts(L, outdir, analysis):
	Scripts = {
		'Taxonomy': Bin / '1.Taxonomy/Taxonomy.py', 
		'ResistanceGene': Bin / '2.ResistanceGene/ResistanceGene.py',
		'VFDB': Bin / '3.VFDB/VFDB.py',
	}
	Workdir = {
		'Taxonomy': outdir / '01.Taxonomy',
		'ResistanceGene': outdir / '02.ResistanceGene',
		'VFDB': outdir / '03.VFDB',
	}
	shell = outdir / 'shell'
	createPath(outdir, shell)
	L2 = ''
	with open(shell / 'Annotation.sh', 'w') as sh:
		for step in analysis.split(','):
			if re.findall('=', step):
				module, version = step.split('=')
				sh.write(f'{cfg.python} {Scripts[module]} --list {L} --db {version} --type {Type} {mode} --outdir {Workdir[module]}\n')
				if module == 'Taxonomy':
					L2 = str(Workdir[module]) + '/sample.list'
			else:
				module = step
				sh.write(f'{cfg.python} {Scripts[module]} --list {L} --type {Type} {mode} --outdir {Workdir[module]}\n')
	if L2:
		return str(shell) + '/Annotation.sh', L2
	else:
		return str(shell) + '/Annotation.sh', L

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
	if mode.upper() == 'SLURM':
		mode = f"--run_mode {mode}"
		if args['partition']:
			mode += f" --partition {args['partition']}"
	else:
		mode = f"--run_mode {mode}"

	cmd, list = generateScripts(infile, outdir, analysis)
	printInfo(os.path.abspath(__file__), 1)
	os.system(f'{cfg.perl} {cfg.watchDog} --num_paral 3 {cmd}')
	os.system(f'ln -s -f {list} {outdir}/sample.list')
	printInfo(os.path.abspath(__file__), 2)
