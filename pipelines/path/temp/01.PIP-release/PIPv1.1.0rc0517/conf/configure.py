#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'qiuwancen'
__email__ = '972538446@qq.com'
__date__ = '2022-04-07'
__version__ = 'V1.0'

import os
import re
import subprocess
from configparser import ConfigParser
Bin = os.path.dirname(os.path.abspath(__file__)) + '/../bin'

config = os.path.abspath(__file__).rstrip('py') + 'ini'
cfg = ConfigParser()
cfg.read(config)
MaxMem = cfg.get('resource', 'MaxMem')
MaxMem = re.sub(r'G', '', MaxMem)
MaxThread = cfg.get('resource', 'MaxThread')
MergeMem = cfg.get('module', 'MergeMem')
MergeMem = re.sub(r'M', '', MergeMem)
FastpMem = cfg.get('module', 'FastpMem')
FastpMem = re.sub(r'G', '', FastpMem)
HostMem = cfg.get('module', 'HostMem')
HostMem = re.sub(r'G', '', HostMem)
HostMemSNAP = cfg.get('module', 'HostMemSNAP')
HostMemSNAP = re.sub(r'G', '', HostMemSNAP)
TaxMem = cfg.get('module', 'TaxMem')
TaxMem = re.sub(r'G', '', TaxMem)
LOCALnum = cfg.get('resource', 'LOCALnum')
SLURMnum = cfg.get('resource', 'SLURMnum')
SGEnum = cfg.get('resource', 'SGEnum')

python = cfg.get('software', 'python')
perl = cfg.get('software', 'perl')
Rscript=cfg.get('software', 'Rscript')
sbatch = cfg.get('software', 'sbatch')
watchDog = cfg.get('software', 'watchDog')
SGE = cfg.get('software', 'SGE')
fastp = cfg.get('software', 'fastp')
sortmerna = cfg.get('software', 'sortmerna')
kraken = cfg.get('software', 'kraken')
krakenMR = cfg.get('software', 'krakenMR')
krakenR = cfg.get('software', 'krakenR')
blastn = cfg.get('software', 'blastn')
bwa2 = cfg.get('software', 'bwa2')
snap = cfg.get('software', 'snap')
bamdst = cfg.get('software', 'bamdst')
samtools = cfg.get('software', 'samtools')
bedtools = cfg.get('software', 'bedtools')
freebayesP = cfg.get('software', 'freebayesP')
pigz = cfg.get('software', 'pigz')
bwa = cfg.get('software', 'bwa')
resfinder = cfg.get('software', 'resfinder')
KMA = cfg.get('software', 'kma')
txt2excel = cfg.get('software', 'txt2excel')
hostDB = cfg.get('database', 'hostDB')
KrakenDB = cfg.get('database', 'KrakenDB')
MEGARes = cfg.get('database', 'MEGARes')
rRNADB = cfg.get('database', 'rRNADB')
resfinderDB = cfg.get('database', 'resfinderDB')
VFDB = cfg.get('database', 'VFDB')
BLASTDB = cfg.get('database', 'BLASTDB')
seqTaxID2Tax = cfg.get('database', 'seqTaxID2Tax')
taxonomy=cfg.get('database', 'taxonomy')

def Resource(mode):
	HostNum = 1
	if mode.upper() == 'SLURM':
		HostNum = int(SLURMnum)
	elif mode.upper() == 'SGE':
		HostNum = int(SGEnum)
	else:
		HostNum = int(LOCALnum)
	MaxThread2 = int(MaxThread) * int(HostNum)
	MaxMem2 = int(MaxMem) * int(HostNum)
	return MaxThread2, MaxMem2

def createPath(*paths):
	for path in paths:
		if not os.path.exists(path):
			os.makedirs(path)

def printInfo(script, flag):
	import time
	if int(flag) == 1:
		print('[%(Time)s] INFO: Start program %(script)s' % \
			{'Time': time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), 'script': script})
	else:
		print('[%(Time)s] INFO: Finish program %(script)s' % \
			{'Time': time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), 'script': script})
	
def checkType(Type, l):
	if (Type == 'SE' and len(l) != 4) or (Type == 'PE' and len(l) != 5):
		print("Error: there is a conflict between the input file format and the sequencing sequence type.")
		exit(0)

def checkFile(infile):
	infile = os.path.abspath(infile)
	if not os.path.isfile(infile):
		print('%(infile)s is not exist!' % {'infile': infile})
		exit(0)

def parseVersion(path):
	with open(path + '/Versions.log', 'w') as out:
		out.write('[Pipeline]\n')
		main_script = os.path.join(Bin, 'PIP')
		_, ver = subprocess.getstatusoutput(f'{main_script} --version')
		out.write(ver + '\n\n')

		out.write('[Software]\n')
		status = ''
		output = ''
		out.write('\n')
		# for root, dirs, files in os.walk(Bin):
		# 	for f in files:
		# 		filePath = os.path.abspath(root + '/' + f)
		# 		(status, output) = subprocess.getstatusoutput(f'{filePath} --version')
		# 		out.write(output + '\n')

		for i in cfg.items('software'):
			if i[0] == 'perl':
				(status, output) = subprocess.getstatusoutput("{} --version |grep version |sed -e 's/).*//g' -e 's/.*(//g' |awk '{}print \"Perl \"$0{}'".format(i[1], '{', '}'))
				out.write(output + '\n')
			elif i[0] == 'fastp':
				(status, output) = subprocess.getstatusoutput("{} 2>&1 |grep version |awk '{}print \"fastp \"$0{}'".format(i[1], '{', '}'))
				out.write(output + '\n')
			elif i[0] == ['python', 'kraken', 'bedtools', 'pigz', 'sbatch', 'watchDog']:
				(status, output) = subprocess.getstatusoutput(f"{i[1]} --version |head -n1")
				out.write(output + '\n')
			elif i[0] == 'bwa2':
				(status, output) = subprocess.getstatusoutput("{} version |awk '{}print \"bwa-mem2 \"$1{}'".format(i[1], '{', '}'))
				out.write(output + '\n')
			elif i[0] == 'bamdst':
				(status, output) = subprocess.getstatusoutput(f"{i[1]} --help |grep version")
				out.write(output + '\n')
			elif i[0] == 'samtools':
				(status, output) = subprocess.getstatusoutput("{} 2>&1 |grep Version |awk '{}print \"SAMtools \"$0{}'".format(i[1], '{', '}'))
				out.write(output + '\n')
			elif i[0] == 'freebayesP':
				(status, output) = subprocess.getstatusoutput(f"echo {i[1]} |sed 's/-.*//g'")
				(status, output) = subprocess.getstatusoutput("{} |grep version |awk '{}print \"FreeBayes \"$0{}'".format(output, '{', '}'))
				out.write(output + '\n')
			elif i[0] == 'bwa':
				(status, output) = subprocess.getstatusoutput("{} 2>&1 |grep Version |awk '{}print \"BWA \"$0{}'".format(i[1], '{', '}'))
				out.write(output + '\n')
		out.write('\n')

# if __name__ == '__main__':
# 	 load_configure()
