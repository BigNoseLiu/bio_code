#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = 'qiuwancen'
__email__ = '972538446@qq.com'
__date__ = '2021-12-04'
__version__ = 'V1.0'

import os
import sys
import math
import pandas as pd
from sys import argv
from pathlib import Path
Bin = Path(__file__).resolve().parent
sys.path.append(os.path.join(Bin, '../../conf'))
import configure as cfg
from configure import createPath
from plotly import offline
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from argparse import ArgumentParser

def readArg(args):
	parser = ArgumentParser(description='{} {} ({}) {}'.format(\
		os.path.basename(os.path.abspath(__file__)), __version__, __date__, __email__))
	args = parser.add_argument
	args('--qc', type=str, action='store', required=True, \
		help='QC table.')
	args('--profile', type=str, action='store', required=True, \
		help='Species profile.')
	args('--outdir', type=str, action='store', default='./', \
		help='Output Path. [%(default)s]')
	args('--version', action='version', version='{} {} ({})'.format(\
		os.path.basename(os.path.abspath(__file__)), __version__, __date__), \
		help='Print the current version and exit.')
	return vars(parser.parse_args())

if __name__ == '__main__':
	args = readArg(argv)
	qc = os.path.abspath(args['qc'])
	profile = os.path.abspath(args['profile'])
	outdir = os.path.abspath(args['outdir'])
	createPath(outdir)

	# common
	#Target = ['Raw Bases Q20', 'Raw Bases Q30', 'Clean Bases Q20', 'Clean Bases Q30', 'Duplication Rate', 'Host Rate', 'Classified_Rate'] #hy
	Target = ['Raw Bases Q20', 'Raw Bases Q30', 'Duplication Rate', 'Host Rate', 'Classified_Rate'] #hy
	Types = ['Eukaryota:Protozoa', 'Viruses', 'Bacteria', 'Eukaryota:Parasite', 'Eukaryota:Fungi']
	Filter = {
		'Eukaryota:Protozoa': 24,
		'Viruses': 18,
		'Bacteria': 24,
		'Eukaryota:Parasite': 24,
		'Eukaryota:Fungi': 24,
	}
	plot_cols = {
		'Eukaryota:Protozoa': 8,
		'Viruses': 6,
		'Bacteria': 8,
		'Eukaryota:Parasite': 8,
		'Eukaryota:Fungi': 8,
	}

	# QC
	df = pd.read_table(qc, header=0, index_col='Sample', delimiter="\t")
	X = list(df.columns)
	Rows = list(df.index)
	fig = go.Figure()
	for i in Target:
		table = [float(k.replace('%', '')) for k in df.iloc[Rows.index(i)]]
		fig.add_trace(go.Scatter(x=X, y=table, name=i))
	offline.plot(fig, filename=f'{outdir}/QC.html',image_width=100, image_height=100, auto_play=False, auto_open=False)

	# Profile
	df = pd.read_table(profile, header=0, index_col='ScientificName', delimiter="\t")
	exList = list(df.Kingdom)
	exList.remove('-')
	for a in Types:
		A = [x for x in exList if x == a]
		df2 = df[df.Kingdom.isin(A)]
		Col = list(df2.columns)
		for b in ['Kingdom', 'ChineseName', 'Comment']:
			Col.remove(b)
		df2 = df2[Col]
		rows = df2.index
		# plotly setup
		drawNum = 72
		if len(rows) < drawNum:
			drawNum = len(rows)
		figNum = math.ceil(drawNum / Filter[a])
		x = 0
		for n in range(1, figNum+1):
			count = 0
			plot_rows = math.ceil(Filter[a] / plot_cols[a])
			fig = make_subplots(rows=plot_rows, cols=plot_cols[a])
			# add traces
			for i in range(1, plot_rows + 1):
				for j in range(1, plot_cols[a] + 1):
					if x < drawNum and count < Filter[a]:
						ha = dict(df2.iloc[x]) #hovertemplate
						fig.add_trace(go.Box(y=[int(k) for k in df2.iloc[x]], name = rows[x], text=Col, notchwidth=0, pointpos=0), row=i, col=j)
						count += 1
						x += 1
			fig.update_traces(boxpoints='all', jitter=1)
			offline.plot(fig, filename=f'{outdir}/{a}.{n}.Boxplot.html',image_width=800, image_height=6000, auto_play=False, auto_open=False)
			os.system(f'cat {outdir}/{a}.{n}.Boxplot.html >> {outdir}/{a}.Boxplot.html')
		os.system(f'rm {outdir}/{a}.*.Boxplot.html')
