#!/bin/sh
#$ -S /bin/sh
perl ../groups_job_list.pl NGB_DB -mem
perl ../qsub-sge.pl --convert no --resource vf=1G --pro ngb_db --queue ngb.q  --jobprefix step1 step1.sh
