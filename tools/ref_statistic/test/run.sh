#!/bin/sh
#$ -S /bin/sh
perl ../ref_fa.stat.pl -ref /ifs1/pub/database/hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa/chrUn_gl000226.fa -out stat.xls
perl ../ref_fa.stat.pl -ref /ifs1/pub/database/hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa/chr1.fa -out stat2.xls
