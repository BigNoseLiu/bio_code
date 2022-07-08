#!/bin/sh
#$ -S /bin/sh
jobNum=1
qstat -u $1|sed '1,2d'|head -$jobNum|awk '{print $1}'|while read first_ID
do
echo "Fisrt job ID: "$first_ID;
qstat -j $first_ID|grep 'usage'
done
