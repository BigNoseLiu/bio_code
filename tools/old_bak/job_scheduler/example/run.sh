#!/bin/sh
#$ -S /bin/sh
mkdir -p log
nohup perl /ifs1/NGB_DB/liumingming/XArchive4Gene/tools/job_scheduler/bin/job_scheduler.pl -s test.schedule -l test.log >run.nohup 2>run.error &
