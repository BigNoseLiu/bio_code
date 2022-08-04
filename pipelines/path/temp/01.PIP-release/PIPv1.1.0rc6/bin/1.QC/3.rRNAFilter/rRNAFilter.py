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
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../conf'))
import configure as cfg
from configure import createPath, printInfo, checkType, checkFile, Resource

def ReadArg(args):
	parser = ArgumentParser(description = '{} {} ({}) {}'.format( \
		os.path.basename(os.path.abspath(__file__)), __version__, __date__, __email__))
	args = parser.add_argument
	args('--list', type=str, action='store', required=True, \
		help='List of Fastq, format: "sample\\tDNA|RNA\\tread1[\\tread2]\\tread_length".')
	args('--type', type=str, action='store', default='SE', choices=['SE', 'PE'], \
		help='Input sequence type. [%(default)s]')
	args('--thread', type=int, action='store', default=4, \
		help='Thread. [%(default)s]')
	args('--run_mode', type=str, action='store', default='local', choices=['local', 'Slurm'], \
		help='Set the run mode, local, Slurm or SGE the option --partition must be set if you set run mode to Slurm or SGE. [%(default)s]')
	args('--partition', type=str, action='store', choices=['ngs-node1', 'ngs-node2', 'ngs-node3'], \
		help='Partition requested, it need to be set if you set --run_mode to Slurm or SGE.')
	args('--outdir', type=str, action='store', default='.', 
		help='Output path. [%(default)s]')
	args('--version', action='version', version='{} {} ({})'.format( \
		os.path.basename(os.path.abspath(__file__)), __version__, __date__), \
		help='Print the current version and exit.')
	return vars(parser.parse_args())

def generateScripts(L, outdir):
	path = Path(outdir)
	shell = path / 'shell'
	createPath(outdir, shell)

	SH = shell / 'rRNAFilter.sh'
	sList = path / 'sample.list'
	line = 0
	with open(L) as inf, open(sList, 'w') as sL, open(SH, 'w') as sh:
		for l in inf.readlines():
			line += 1
			l = l.strip().split('\t')
			checkType(Type, l)

			samplePath = path / l[0] / l[1]
			createPath(samplePath)

			if Type == 'SE':
				sh.write('rm -fr %(samplePath)s/kvdb\n'
					'%(sortmerna)s --threads %(thread)s --ref %(rRNADB)s/smr_db.fasta --reads %(rawRead)s '
					'--workdir %(samplePath)s --idx %(rRNADB)s/idx --fastx --aligned %(samplePath)s/%(sample)s.rRNA '
					'--other %(samplePath)s/%(sample)s.no_rRNA --num_alignments 1\n'
					'%(pigz)s -p %(thread)s %(samplePath)s/%(sample)s.rRNA.fq\n'
					'%(pigz)s -p %(thread)s %(samplePath)s/%(sample)s.no_rRNA.fq\n'% \
					{'sortmerna': cfg.sortmerna, 'rRNADB': cfg.rRNADB, 'thread': args['thread'], 'rawRead': l[2], 
					'samplePath': samplePath, 'sample': l[0], 'pigz': cfg.pigz})
				sL.write(f'{l[0]}\t{l[1]}\t{samplePath}/{l[0]}.no_rRNA.fq.gz\t{l[-1]}\n')
			else:
				sh.write('rm -fr %(samplePath)s/kvdb\n'
					'%(sortmerna)s --threads %(thread)s --ref %(rRNADB)s/smr_db.fasta --reads %(rawRead1)s --reads %(rawRead2)s '
					'--workdir %(samplePath)s --idx %(rRNADB)s/idx --fastx --aligned %(samplePath)s/%(sample)s.rRNA '
					'--other %(samplePath)s/%(sample)s.no_rRNA --num_alignments 1 --paired_in --out2\n'
					'%(pigz)s -p %(thread)s %(samplePath)s/%(sample)s.rRNA_fwd.fq\n'
					'%(pigz)s -p %(thread)s %(samplePath)s/%(sample)s.rRNA_rev.fq\n'
					'%(pigz)s -p %(thread)s %(samplePath)s/%(sample)s.no_rRNA_fwd.fq\n'
					'%(pigz)s -p %(thread)s %(samplePath)s/%(sample)s.no_rRNA_rev.fq\n' % \
					{'sortmerna': cfg.sortmerna, 'rRNADB': cfg.rRNADB, 'thread': args['thread'], 'rawRead1': l[2], 'rawRead2': l[3], 
					'samplePath': samplePath, 'sample': l[0], 'pigz': cfg.pigz})
				sL.write(f'{l[0]}\t{l[1]}\t{samplePath}/{l[0]}.no_rRNA_fwd.fq.gz\t{samplePath}/{l[0]}.no_rRNA_rev.fq.gz\t{l[-1]}\n')
	return SH, line


if __name__ == '__main__':
	args = ReadArg(argv)
	infile = os.path.abspath(args['list'])
	outdir = os.path.abspath(args['outdir'])
	Type = args['type']
	mode = args['run_mode']

	checkFile(infile)

#	if mode.upper() == 'SLURM' and not args['partition']:
#		print('Error: the option --partition must be set when you set --run_mode to Slurm')
#		exit(0)

	num_paral = len(open(f'{infile}').readlines())
	MaxThread, MaxMem = Resource(mode)

	SH, lineNum = generateScripts(infile, outdir)
	printInfo(os.path.abspath(__file__), 1)
	if mode.upper() == 'SLURM':
		if Type == 'SE':
			os.system(f'{cfg.perl} {cfg.sbatch} --maxjob 105 --lines 4 --cpu {args["thread"]} --jobprefix rRNAFilter {SH}')
		else:
			os.system(f'{cfg.perl} {cfg.sbatch} --maxjob 105 --lines 6 --cpu {args["thread"]} --jobprefix rRNAFilter {SH}')
	else:
		if Type == 'SE':
			os.system(f'{cfg.perl} {cfg.watchDog} --mem 8G --num_paral 35 --num_line 4 {SH}')
		else:
			os.system(f'{cfg.perl} {cfg.watchDog} --mem 8G --num_paral 35 --num_line 6 {SH}')
	printInfo(os.path.abspath(__file__), 2)
