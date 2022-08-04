#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'qiuwancen'
__email__ = '972538446@qq.com'
__date__ = '2022-04-14'
__version__ = 'V1.1'

import os
import sys
from sys import argv
from pathlib import Path
from argparse import ArgumentParser

Bin = os.path.dirname(os.path.abspath(__file__))

sys.path.append(os.path.join(Bin, '../../../conf'))
import configure as cfg
from configure import createPath, checkFile, checkType, printInfo, Resource

###################### scripts ######################
krakenFilter = os.path.join(Bin, 'krakenFilter.py')
TaxSplitGenus = os.path.join(Bin, 'TaxSplitGenus.pl')
getSpeciesReads = os.path.join(Bin, 'getSpeciesReads.pl')
reportUpdate = os.path.join(Bin, 'reportUpdate.pl')
TaxonomyStat = os.path.join(Bin, 'TaxonomyStat.py')
Result2BWA = os.path.join(Bin, 'Result2BWA.pl')
Ufastq = os.path.join(Bin, 'Ufastq.pl')
KoutScore = os.path.join(Bin, 'KoutScore.py')
Mark2Fa = os.path.join(Bin, 'Mark2Fa.pl')
bwaAlign = os.path.join(Bin, 'bwaAlign.py')
coverageStat = os.path.join(Bin, 'coverageStat.py')
coveragePlot = os.path.join(Bin, 'coverage2.R')
#####################################################

def readArg(args):
	parser = ArgumentParser(description = \
		f'{os.path.basename(os.path.abspath(__file__))} {__version__} ({__date__}) {__email__}')
	args = parser.add_argument
	args('--list', type=str, action='store', required=True, \
		help='List of Fastq, format: "sample\\tDNA|RNA\\tread1[\\tread2]\\tsample type", sample type contain: RT,CSF,Plasma,[Others].')
	args('--db', type=str, action='store', default='PIDB202101_Middle_addHPV18_taxidCorr', \
		help='Database version of PIP. [%(default)s]')
	args('--thread', type=int, action='store', default=48, \
		help='Thread. [%(default)s]')
	args('--confidence', type=float, action='store', default=0.0, \
		help='Confidence score threshold, must be in [0, 1]. [%(default)s]')
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

	dbVersion = args['db']

	SH = shell / 'Kraken.sh'
	SH2 = shell / 'Coverage.sh'
	sList = outdir / 'sample.list'
	line = 0
	with open(L) as inf, open(sList, 'w') as sL, open(SH, 'w') as sh, open(SH2, 'w') as sh2:
		for fullLine in inf.readlines():
			line += 1
			fullLine = fullLine.strip()
			l = fullLine.split('\t')
			checkType(Type, l)

			samplePath = outdir / l[0] / l[1]
			createPath(samplePath)

			if l[2].endswith('gz'):
				gzip = '--gzip-compressed '
			else:
				gzip = ''

			if Type == 'SE':
				sh.write('%(perl)s %(kraken)s --db %(KrakenDB)s/%(dbVersion)s --confidence %(confidence)s --threads 48 %(gzip)s'
					'--report %(samplePath)s/%(sample)s.raw.kreport --classified-out %(samplePath)s/%(sample)s.kout.C.fastq '
					'--output %(samplePath)s/%(sample)s.raw.kout %(fastq)s\n' % \
					{'kraken': cfg.kraken, 'KrakenDB': cfg.KrakenDB, 'dbVersion': args['db'], 'confidence': args['confidence'], 'gzip': gzip, \
					'samplePath': samplePath, 'sample': l[0], 'fastq': l[2], 'krakenMR': cfg.krakenMR, 'krakenR': cfg.krakenR, 'perl': cfg.perl})
			else:
				sh.write('%(perl)s %(kraken)s --db %(KrakenDB)s/%(dbVersion)s --confidence %(confidence)s --threads 48 %(gzip)s'
					'--report %(samplePath)s/%(sample)s.raw.kreport --classified-out %(samplePath)s/%(sample)s.kout.C.fastq#.fq '
					'--output %(samplePath)s/%(sample)s.raw.kout --paired %(fastq1)s %(fastq2)s\n' % \
					{'kraken': cfg.kraken, 'KrakenDB': cfg.KrakenDB, 'dbVersion': args['db'], 'confidence': args['confidence'], 'gzip': gzip, \
					'samplePath': samplePath, 'sample': l[0], 'fastq1': l[2], 'fastq2': l[3], 'perl': cfg.perl})
			sh.write('%(python)s %(krakenFilter)s --ko %(samplePath)s/%(sample)s.raw.kout --db %(dbVersion)s --fc unique --out %(samplePath)s/%(sample)s.unique.kout\n'
					'%(perl)s %(krakenMR)s --db %(KrakenDB)s/%(dbVersion)s %(samplePath)s/%(sample)s.unique.kout > %(samplePath)s/%(sample)s.unique.mpa.kreport &\n'
					'%(perl)s %(krakenR)s --db %(KrakenDB)s/%(dbVersion)s %(samplePath)s/%(sample)s.unique.kout > %(samplePath)s/%(sample)s.unique.kreport &\n'
					'%(python)s %(KoutScore)s --kout %(samplePath)s/%(sample)s.raw.kout --top 5 --db %(KrakenDB)s/%(dbVersion)s --out %(samplePath)s/%(sample)s.raw.kout.MarkNum &\n'
					'wait\n' % \
					{'krakenFilter': krakenFilter, 'samplePath': samplePath, 'sample': l[0], 'KrakenDB': cfg.KrakenDB, 'KoutScore': KoutScore, \
					'dbVersion': args['db'], 'Rlen': l[-1], 'krakenMR': cfg.krakenMR, 'krakenR': cfg.krakenR, 'python': cfg.python, 'perl': cfg.perl})
			if Type == 'SE':
				sh.write('%(perl)s %(TaxSplitGenus)s --report %(samplePath)s/%(sample)s.unique.kreport --db %(KrakenDB)s/%(dbVersion)s --fq %(fastq)s '
					'--inMark %(samplePath)s/%(sample)s.raw.kout.MarkNum --krakenout %(samplePath)s/%(sample)s.unique.kout --outdir %(samplePath)s\n' % \
#				sh.write('%(perl)s %(Mark2Fa)s --inMark %(samplePath)s/%(sample)s.raw.kout.MarkNum --db 20210103 '
#					'--fq %(samplePath)s/%(sample)s.kout.C.fastq --outdir %(samplePath)s\n' % \
					{'TaxSplitGenus': TaxSplitGenus, 'fastq': l[2], 'samplePath': samplePath, 'sample': l[0], 'KrakenDB': cfg.KrakenDB, \
					'dbVersion': args['db'], 'perl': cfg.perl, 'Mark2Fa': Mark2Fa})
			else:
				sh.write('%(perl)s %(TaxSplitGenus)s --report %(samplePath)s/%(sample)s.unique.kreport --db %(KrakenDB)s/%(dbVersion)s --fq %(fastq1)s,%(fastq2)s ' 
					'--inMark %(samplePath)s/%(sample)s.raw.kout.MarkNum --krakenout %(samplePath)s/%(sample)s.unique.kout --outdir %(samplePath)s\n' % \
#				sh.write('%(perl)s %(Mark2Fa)s --inMark %(samplePath)s/%(sample)s.raw.kout.MarkNum --db 20210103 '
#					'--fq %(samplePath)s/%(sample)s.kout.C.fastq_1.fq,%(samplePath)s/%(sample)s.kout.C.fastq_2.fq --outdir %(samplePath)s\n' % \
					{'TaxSplitGenus': TaxSplitGenus, 'fastq1': l[2], 'fastq2': l[3], 'samplePath': samplePath, 'sample': l[0],
					'KrakenDB': cfg.KrakenDB, 'dbVersion': args['db'], 'perl': cfg.perl, 'Mark2Fa': Mark2Fa})
			sh.write('%(perl)s %(reportUpdate)s --raw %(samplePath)s/%(sample)s.raw.kreport --update %(samplePath)s/%(sample)s.unique.kreport\n'
			'%(python)s %(TaxonomyStat)s --report %(samplePath)s/%(sample)s.unique.kreport --mpa %(samplePath)s/%(sample)s.unique.mpa.kreport --type %(type)s '
			'--db %(KrakenDB)s/%(dbVersion)s --name %(sample)s_%(moleculeType)s --validate %(samplePath)s/validated.out '
			'--highconf %(samplePath)s/Species-HighConfidenceKmerScore.txt --outdir %(samplePath)s\n'
			f'{cfg.python} {bwaAlign} -i {samplePath} --db {cfg.KrakenDB}/{dbVersion}/BWA.v5 -o {samplePath}\n'
			'%(perl)s %(Result2BWA)s --cmd %(samplePath)s/BWAcmd.list --pathoList %(samplePath)s/Pathogeny_Detail.txt --outdir %(samplePath)s\n' % \
			{'samplePath': samplePath, 'sample': l[0], 'KrakenDB': cfg.KrakenDB, 'dbVersion': args['db'], 'type': l[-1], \
                        'reportUpdate': reportUpdate, 'TaxonomyStat': TaxonomyStat, 'moleculeType': l[1], 'perl': cfg.perl, 'Result2BWA': Result2BWA, 'python': cfg.python})
			sh2.write(f'{cfg.python} {coverageStat} -s {coveragePlot} -i {samplePath}/BwaAlign/ --db {cfg.KrakenDB}/{dbVersion}/BWA.v5 -p 10 -o {samplePath}/CovPlot\n')
			sL.write(f'{fullLine}\t{samplePath}/{l[0]}.unique.kout\t{samplePath}/{l[0]}.unique.kreport\t')
			sL.write(f'{samplePath}/{l[0]}.unique.mpa.kreport\t{samplePath}/taxid.list\t{samplePath}/Total_Detail.txt\t{samplePath}/Pathogeny_Detail.txt\n')
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

	num_paral = len(open(f'{infile}').readlines())
	MaxThread, MaxMem = Resource(mode)
	TaxMem = int(os.path.getsize(f'{cfg.KrakenDB}/{args["db"]}/hash.k2d') / 1024 / 1024 / 1024) + 1
	num_paral = int(MaxMem / TaxMem)
	thread = int(MaxThread / num_paral) + 1
	if thread > 24:
		thread = 24

	cmd, lineNum = generateScripts(infile, outdir)
	#exit(0)
	printInfo(os.path.abspath(__file__), 1)
	if mode.upper() == 'SLURM':
		os.system(f'{cfg.perl} {cfg.sbatch} --maxjob {num_paral} --cpu {thread} --lines 11 --jobprefix Taxonomy {cmd}')
		os.system(f'cat {outdir}/*/*/BWA.sh > {outdir}/shell/BWA.sh')
		os.system(f'{cfg.perl} {cfg.sbatch} --maxjob {MaxThread} --cpu 1 --jobprefix BWA {outdir}/shell/BWA.sh')
		os.system(f'{cfg.perl} {cfg.sbatch} --maxjob {MaxThread} --cpu 1 --jobprefix Coverage {outdir}/shell/Coverage.sh')
	elif mode.upper() == 'SGE':
		os.system(f'{cfg.perl} {cfg.SGE} --resource vf={TaxMem}g,p={thread} --lines 11 --queue all.q --maxjob {num_paral} --jobprefix Taxonomy {cmd}')
		os.system(f'cat {outdir}/*/*/BWA.sh > {outdir}/shell/BWA.sh')
		os.system(f'{cfg.perl} {cfg.SGE} --resource vf=2g,p=1 --queue all.q --maxjob {MaxThread} --jobprefix BWA {outdir}/shell/BWA.sh')
		os.system(f'{cfg.perl} {cfg.SGE} --resource vf=2g,p=1 --queue all.q --maxjob {MaxThread} --jobprefix Coverage {outdir}/shell/Coverage.sh')
	else:
		os.system(f'{cfg.perl} {cfg.watchDog} --mem {TaxMem}G --num_paral {num_paral} --num_line 11 {cmd}')
		os.system(f'cat {outdir}/*/*/BWA.sh > {outdir}/shell/BWA.sh')
		os.system(f'{cfg.perl} {cfg.watchDog} --mem 2G --num_paral {MaxThread} {outdir}/shell/BWA.sh')
		os.system(f'{cfg.perl} {cfg.watchDog} --mem 2G --num_paral {MaxThread} {outdir}/shell/Coverage.sh')
		pass
	printInfo(os.path.abspath(__file__), 2)
