#!/bin/sh
#$ -S /bin/sh
mkdir shell log original_stat

ls /ifs1/pub/database/hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa/|while read L
do
	echo "#!/bin/sh
#$ -S /bin/sh
perl /ifs1/NGB_DB/liumingming/XArchive4Gene/tools/ref_statistic/ref_fa.stat.pl -ref /ifs1/pub/database/hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa/$L -out original_stat/$L.original_stat
">shell/stat_$L.sh
done

ls shell/*|while read sh
do
	#qsub -cwd -l vf=200m -q ngb.q -P ngb_db -o log/ -e log/ $sh
done
