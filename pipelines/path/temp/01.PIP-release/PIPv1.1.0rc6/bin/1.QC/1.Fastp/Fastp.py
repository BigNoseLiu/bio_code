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
	args('--thread', type=int, action='store', default=8, \
		help='Thread. [%(default)s]')
	args('--adaptor3', type=str, action='store', default='AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA', \
		help='The adapter for read1. [%(default)s]')
	args('--adaptor5', type=str, action='store', default='AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG', \
		help='The adapter for read2 (PE data only). [%(default)s]')
	args('--run_mode', type=str, action='store', default='local', choices=['local', 'Slurm', 'SGE'], \
		help='Set the run mode, local, Slurm or SGE, the option --partition must be set if you set run mode to Slurm or SGE. [%(default)s]')
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

	SH = shell / 'Fastp.sh'
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
				sh.write('%(fastp)s -y -q 20 -u 50 -n 5 --length_required 50 --adapter_sequence %(adaptor3)s -c --thread %(thread)s -i %(rawRead)s -o %(samplePath)s/%(sample)s.fq.gz '
					'-j %(samplePath)s/%(sample)s.json -h %(samplePath)s/%(sample)s.html 2>%(samplePath)s/%(sample)s.log\n' % \
					{'fastp': cfg.fastp, 'adaptor3': args['adaptor3'], 'thread': thread, 'rawRead': l[2], 'samplePath': samplePath, 'sample': l[0], 'libType': l[1]})
				sL.write(f'{l[0]}\t{l[1]}\t{samplePath}/{l[0]}.fq.gz\t{l[-1]}\n')
			else:
				sh.write('%(fastp)s -y -q 20 -u 50 -n 5 --length_required 50 --adapter_sequence %(adaptor3)s --adapter_sequence_r2 %(adaptor5)s -c --thread %(thread)s -i %(rawRead1)s -I %(rawRead2)s '
					'-o %(samplePath)s/%(sample)s_1.fq.gz -O %(samplePath)s/%(sample)s_2.fq.gz '
					'-j %(samplePath)s/%(sample)s.json -h %(samplePath)s/%(sample)s.html 2>%(samplePath)s/%(sample)s.log\n' % \
					{'fastp': cfg.fastp, 'adaptor3': args['adaptor3'], 'adaptor5': args['adaptor5'], 'thread': thread, 'rawRead1': l[2], 'rawRead2': l[3], 'samplePath': samplePath, 'sample': l[0], 'libType': l[1]})
				sL.write(f'{l[0]}\t{l[1]}\t{samplePath}/{l[0]}_1.fq.gz\t{samplePath}/{l[0]}_2.fq.gz\t{l[-1]}\n')
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
	if int(num_paral) > MaxThread:
		num_paral = MaxThread
	if num_paral > int(MaxMem / int(cfg.FastpMem)):
		num_paral = int(MaxMem / int(cfg.FastpMem))
	jobThread = int(int(MaxThread) / num_paral)
	thread = jobThread + 1
	if jobThread > 8:
		jobThread = 8
		thread = 8
	
	SH, lineNum = generateScripts(infile, outdir)

	printInfo(os.path.abspath(__file__), 1)
	if mode.upper() == 'SLURM':
		os.system(f'{cfg.perl} {cfg.sbatch} --maxjob {num_paral} --cpu {jobThread} --jobprefix Fastp {SH}')
	elif mode.upper() == 'SGE':
		os.system(f'{cfg.perl} {cfg.SGE} --resource vf={cfg.FastpMem}g,p={jobThread} --queue all.q --maxjob {num_paral} --jobprefix Fastp {SH}')
	else:
		os.system(f'{cfg.perl} {cfg.watchDog} --mem {cfg.FastpMem}G --num_paral {num_paral} {SH}')
	printInfo(os.path.abspath(__file__), 2)
