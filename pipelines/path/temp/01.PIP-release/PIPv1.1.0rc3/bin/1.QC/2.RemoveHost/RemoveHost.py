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
	args('--soft', type=str, action='store', default='snap-aligner', choices=['bwa-mem2', 'snap-aligner'], \
		help='Align soft. [%(default)s]')
	args('--thread', type=int, action='store', default=20, \
		help='Thread. [%(default)s]')
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

	SH = shell / 'RemoveHost.sh'
	sList = path / 'sample.list'
	line = 0
	with open(L) as inf, open(sList, 'w') as sL, open(SH, 'w') as sh:
		for l in inf.readlines():
			line += 1
			l = l.strip().split('\t')
			checkType(Type, l)

			samplePath = path / l[0] / l[1]
			createPath(samplePath) 

			if soft == 'bwa-mem2':
				if Type == 'SE':
					sh.write('%(bwa)s mem -t %(thread)s -Y -M %(hostDB)s/human.fna %(Reads)s | '
						'%(samtools)s view -@ %(thread)s -bS -o %(samplePath)s/%(sample)s.bam -\n' % \
						{'bwa': INFO[soft]['soft'], 'thread': thread, 'hostDB': INFO[soft]['db'], 'Reads': l[2],
						'samtools': cfg.samtools, 'samplePath': samplePath, 'sample': l[0]})
					sh.write('%(samtools)s fastq -f 4 -@ %(thread)s %(samplePath)s/%(sample)s.bam |gzip > %(samplePath)s/%(sample)s.fq.gz &\n'
						'%(samtools)s flagstat -@ %(thread)s %(samplePath)s/%(sample)s.bam > %(samplePath)s/%(sample)s.Homo.log &\n' 
						'wait\n' % \
						{'thread': int(thread/2), 'samtools': cfg.samtools, 'samplePath': samplePath, 'sample': l[0]}) 
					sL.write(f'{l[0]}\t{l[1]}\t{samplePath}/{l[0]}.fq.gz\t{l[-1]}\n')
				else:
					sh.write('%(bwa)s mem -t %(thread)s -Y -M %(hostDB)s/human.fna %(Read1)s %(Read2)s | '
						'%(samtools)s view -@ %(thread)s -bS -o %(samplePath)s/%(sample)s.bam -\n' % \
						{'bwa': INFO[soft]['soft'], 'thread': thread, 'hostDB': INFO[soft]['db'], 'Read1': l[2], 
						'Read2': l[3], 'samtools': cfg.samtools, 'samplePath': samplePath, 'sample': l[0]})
					sh.write('%(samtools)s fastq -f 12 -@ %(thread)s -1 %(samplePath)s/%(sample)s_1.fq.gz '
						'-2 %(samplePath)s/%(sample)s_2.fq.gz %(samplePath)s/%(sample)s.bam &\n' 
						'%(samtools)s flagstat -@ %(thread)s %(samplePath)s/%(sample)s.bam > %(samplePath)s/%(sample)s.Homo.log &\n' 
						'wait\n' % \
						{'thread': int(thread/2), 'samtools': cfg.samtools, 'samplePath': samplePath, 'sample': l[0]})
					sL.write(f'{l[0]}\t{l[1]}\t{samplePath}/{l[0]}_1.fq.gz\t{samplePath}/{l[0]}_2.fq.gz\t{l[-1]}\n')
			elif soft == 'snap-aligner':
				if Type == 'SE':
					sh.write('%(snap)s single %(db)s -t %(thread)s -compressedFastq %(Reads)s -o -bam %(samplePath)s/%(sample)s.bam '
						'> %(samplePath)s/%(sample)s.SNAP.log\n' % \
						{'snap': INFO[soft]['soft'], 'thread': thread, 'db': INFO[soft]['db'], 'Reads': l[2],
						'samplePath': samplePath, 'sample': l[0]})
					sh.write('%(samtools)s fastq -f 4 -@ %(thread)s %(samplePath)s/%(sample)s.bam |gzip > %(samplePath)s/%(sample)s.fq.gz &\n'
						'%(samtools)s flagstat -@ %(thread)s %(samplePath)s/%(sample)s.bam > %(samplePath)s/%(sample)s.Homo.log &\n'
						'wait\n' % \
						{'thread': int(thread/2), 'samtools': cfg.samtools, 'samplePath': samplePath, 'sample': l[0]})
					sL.write(f'{l[0]}\t{l[1]}\t{samplePath}/{l[0]}.fq.gz\t{l[-1]}\n')
				else:
					sh.write('%(snap)s paired %(db)s -t %(thread)s -compressedFastq %(Read1)s %(Read2)s -o -bam %(samplePath)s/%(sample)s.bam '
						'> %(samplePath)s/%(sample)s.SNAP.log\n' % \
						{'snap': INFO[soft]['soft'], 'thread': thread, 'db': INFO[soft]['db'], 'Read1': l[2],
						'Read2': l[3], 'samplePath': samplePath, 'sample': l[0]})
					sh.write('%(samtools)s fastq -f 12 -@ %(thread)s -1 %(samplePath)s/%(sample)s_1.fq.gz '
						'-2 %(samplePath)s/%(sample)s_2.fq.gz %(samplePath)s/%(sample)s.bam &\n'
						'%(samtools)s flagstat -@ %(thread)s %(samplePath)s/%(sample)s.bam > %(samplePath)s/%(sample)s.Homo.log &\n'
						'wait\n' % \
						{'thread': int(thread/2), 'samtools': cfg.samtools, 'samplePath': samplePath, 'sample': l[0]})
					sL.write(f'{l[0]}\t{l[1]}\t{samplePath}/{l[0]}_1.fq.gz\t{samplePath}/{l[0]}_2.fq.gz\t{l[-1]}\n')
	return SH, line

if __name__ == '__main__':
	args = ReadArg(argv)
	infile = os.path.abspath(args['list'])
	outdir = os.path.abspath(args['outdir'])
	Type = args['type']
	mode = args['run_mode']
	soft = args['soft']

	checkFile(infile)

	INFO = {
		"bwa-mem2": 
			{
				"mem": cfg.HostMem,
				"db": cfg.hostDB,
				"soft": cfg.bwa2,
			},
		"snap-aligner": 
			{
				"mem": cfg.HostMemSNAP,
				"db": f"{cfg.hostDB}/SNAP2",
				"soft": cfg.snap,
			},
	}

	num_paral = len(open(f'{infile}').readlines())
	MaxThread, MaxMem = Resource(mode)
	num_paral = int(int(MaxMem) / int(INFO[soft]['mem']))
	thread = int(MaxThread / num_paral) + 1

	SH, lineNum = generateScripts(infile, outdir)

	printInfo(os.path.abspath(__file__), 1)
	if mode.upper() == 'SLURM':
		os.system(f'{cfg.perl} {cfg.sbatch} --maxjob {num_paral} --cpu {thread} --lines 4 --jobprefix RemoveHost {SH}')
	elif mode.upper() == 'SGE':
		os.system(f'{cfg.perl} {cfg.SGE} --resource vf={INFO[soft]["mem"]}g,p={thread} --queue all.q --lines 4 --maxjob {num_paral} --jobprefix RemoveHost {SH}')
	else:
		os.system(f'{cfg.perl} {cfg.watchDog} --mem {INFO[soft]["mem"]}G --num_paral {num_paral} --num_line 4 {SH}')
	printInfo(os.path.abspath(__file__), 2)
