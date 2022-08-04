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

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../conf'))
import configure as cfg
from configure import createPath, printInfo, checkType, checkFile, Resource

def readArg(args):
	parser = ArgumentParser(description = '{} {} ({}) {}'.format( \
		os.path.basename(os.path.abspath(__file__)), __version__, __date__, __email__))
	args = parser.add_argument
	args('--list', type=str, action='store', required=True, \
		help='List of Fastq, format: "sample\\tDNA|RNA\\tread1[\\tread2]\\tread_length".')
	args('--type', type=str, action='store', default='SE', choices=['SE', 'PE'], \
		help='Input sequence type. [%(default)s]')
	args('--run_mode', type=str, action='store', default='local', choices=['local', 'Slurm', 'SGE'], \
		help='Set the run mode, local, Slurm or SGE, the option --partition must be set if you set run mode to Slurm or SGE. [%(default)s]')
	args('--partition', type=str, action='store', choices=['ngs-node1', 'ngs-node2', 'ngs-node3'], \
		help='partition requested, it need to be set if you set --run_mode to Slurm or SGE.')
	args('--outdir', type=str, action='store', default='.', \
		help='Output path. [%(default)s]')
	args('--version', action='version', version='{} {} ({})'.format( \
		os.path.basename(os.path.abspath(__file__)), __version__, __date__), \
		help='Print the current version and exit.')
	return vars(parser.parse_args())

def generateScripts(L, outdir):
	shell = os.path.join(outdir, 'shell')
	createPath(outdir, shell)
	Tab = {}
	with open(L) as inf:
		for l in inf.readlines():
			l = l.strip().split('\t')
			if not l[0] or l[0].startswith('#'):
				continue

			checkType(Type, l)

			samplePath = os.path.join(outdir, l[0], l[1])
			createPath(samplePath)
			if l[0] not in Tab:
				Tab.update({l[0]: {l[1]: {'read1': [], 'length': int()}}})
				if Type == 'PE':
					Tab[l[0]][l[1]].update({'read2': []})
			if l[1] not in Tab[l[0]]:
				Tab[l[0]].update({l[1]: {'read1': [], 'length': int()}})
				if Type == 'PE':
					Tab[l[0]][l[1]].update({'read2': []})
			Tab[l[0]][l[1]]['read1'].append(l[2])
			Tab[l[0]][l[1]]['length'] = l[-1]
			if Type == 'PE':
				Tab[l[0]][l[1]]['read2'].append(l[3])
	
	with open(shell + '/Merge.sh', 'w') as sh, open(outdir + '/sample.list', 'w') as sL:
		for s in Tab:
			for t in Tab[s]:
				if len(Tab[s][t]['read1']) == 1:
					if Tab[s][t]['read1'][0].endswith('gz'):
						sh.write('ln -s -f %(read1)s %(outdir)s/%(sample)s/%(libType)s/%(sample)s_1.fq.gz\n' % \
							{'read1': Tab[s][t]['read1'][0], 'outdir': outdir, 'sample': s, 'libType': t})
						sL.write('%(sample)s\t%(libType)s\t%(outdir)s/%(sample)s/%(libType)s/%(sample)s_1.fq.gz' % \
							{'outdir': outdir, 'sample': s, 'libType': t})
					else:
						sh.write('ln -s -f %(read1)s %(outdir)s/%(sample)s/%(libType)s/%(sample)s_1.fq\n' % \
							{'read1': Tab[s][t]['read1'][0], 'outdir': outdir, 'sample': s, 'libType': t})
						sL.write('%(sample)s\t%(libType)s\t%(outdir)s/%(sample)s/%(libType)s/%(sample)s_1.fq' % \
							{'outdir': outdir, 'sample': s, 'libType': t})
					if Type == 'PE':
						if Tab[s][t]['read2'][0].endswith('gz'):
							sh.write('ln -s -f %(read2)s %(outdir)s/%(sample)s/%(libType)s/%(sample)s_2.fq.gz\n' % \
								{'read2': Tab[s][t]['read2'][0], 'outdir': outdir, 'sample': s, 'libType': t})
							sL.write('\t%(outdir)s/%(sample)s/%(libType)s/%(sample)s_2.fq.gz' % \
								{'outdir': outdir, 'sample': s, 'libType': t})
						else:
							sh.write('ln -s -f %(read2)s %(outdir)s/%(sample)s/%(libType)s/%(sample)s_2.fq\n' % \
								{'read2': Tab[s][t]['read2'][0], 'outdir': outdir, 'sample': s, 'libType': t})
							sL.write('\t%(outdir)s/%(sample)s/%(libType)s/%(sample)s_2.fq' % \
								{'outdir': outdir, 'sample': s, 'libType': t})
					sL.write('\t%(len)s\n' % {'len': Tab[s][t]['length']})
				else:
					read1 = ' '.join(Tab[s][t]['read1'])
					if Tab[s][t]['read1'][0].endswith('gz'):
						sh.write('cat %(read1)s > %(outdir)s/%(sample)s/%(libType)s/%(sample)s_1.fq.gz\n' % \
							{'read1': read1, 'outdir': outdir, 'sample': s, 'libType': t})
						sL.write('%(sample)s\t%(libType)s\t%(outdir)s/%(sample)s/%(libType)s/%(sample)s_1.fq.gz' % \
							{'outdir': outdir, 'sample': s, 'libType': t})
					else:
						sh.write('cat %(read1)s > %(outdir)s/%(sample)s/%(libType)s/%(sample)s_1.fq\n' % \
							{'read1': read1, 'outdir': outdir, 'sample': s, 'libType': t})
						sL.write('%(sample)s\t%(libType)s\t%(outdir)s/%(sample)s/%(libType)s/%(sample)s_1.fq' % \
							{'outdir': outdir, 'sample': s, 'libType': t})
					if Type == 'PE':
						read2 = ' '.join(Tab[s][t]['read2'])
						if Tab[s][t]['read2'][0].endswith('gz'):
							sh.write('cat %(read2)s > %(outdir)s/%(sample)s/%(libType)s/%(sample)s_2.fq.gz\n' % \
								{'read2': read2, 'outdir': outdir, 'sample': s, 'libType': t})
							sL.write('\t%(outdir)s/%(sample)s/%(libType)s/%(sample)s_2.fq.gz' % \
								{'outdir': outdir, 'sample': s, 'libType': t})
						else:
							sh.write('cat %(read2)s > %(outdir)s/%(sample)s/%(libType)s/%(sample)s_2.fq\n' % \
								{'read2': read2, 'outdir': outdir, 'sample': s, 'libType': t})
							sL.write('\t%(outdir)s/%(sample)s/%(libType)s/%(sample)s_2.fq' % \
								{'outdir': outdir, 'sample': s, 'libType': t})
					sL.write('\t%(len)s\n' % {'len': Tab[s][t]['length']})
	return shell + '/Merge.sh'

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

	num_paral = len(open(f'{infile}').readlines())
	MaxThread, MaxMem = Resource(mode)
	if int(num_paral) > MaxThread:
		num_paral = MaxThread
	if num_paral > int(MaxMem * 1024 / int(cfg.MergeMem)):
		num_paral = int(MaxMem * 1024 / int(cfg.MergeMem))

	cmd = generateScripts(infile, outdir)
	printInfo(os.path.abspath(__file__), 1)
	if mode.upper() == 'SLURM':
		os.system(f'{cfg.perl} {cfg.sbatch} --maxjob {num_paral} --jobprefix Merge {cmd}')
	elif mode.upper() == 'SGE':
		os.system(f'{cfg.perl} {cfg.SGE} --resource vf={cfg.MergeMem}M,p=1 --queue all.q --maxjob {num_paral} --jobprefix Merge {cmd}')
	else:
		os.system(f'{cfg.perl} {cfg.watchDog} --mem {cfg.MergeMem}M --num_paral {num_paral} {cmd}')
	printInfo(os.path.abspath(__file__), 2)
