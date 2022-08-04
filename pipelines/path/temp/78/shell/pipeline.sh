#!/usr/bin/bash
################## Check Tasks
/data/Software/Anaconda/anaconda3/bin/python /data/Pipeline/01.PIP-release/PIPv1.1.0rc7/bin/4.Stat/CheckTasks.py --samples /data/application/pathongen/data/result/V350071671_20220728163235/78/sample.list --workdir /data/application/pathongen/data/result/V350071671_20220728163235/78 &
################## Merge
/data/Software/Anaconda/anaconda3/bin/python /data/Pipeline/01.PIP-release/PIPv1.1.0rc7/bin/0.Merge/Merge.py --list /data/application/pathongen/data/result/V350071671_20220728163235/78/sample.list --type SE --run_mode local --outdir /data/application/pathongen/data/result/V350071671_20220728163235/78/00.Merge
################## QC
/data/Software/Anaconda/anaconda3/bin/python /data/Pipeline/01.PIP-release/PIPv1.1.0rc7/bin/1.QC/QC.py --list /data/application/pathongen/data/result/V350071671_20220728163235/78/00.Merge/sample.list --type SE --run_mode local --outdir /data/application/pathongen/data/result/V350071671_20220728163235/78/01.QC --analysis Fastp,RemoveHost
################## Annotation
/data/Software/Anaconda/anaconda3/bin/python /data/Pipeline/01.PIP-release/PIPv1.1.0rc7/bin/2.Annotation/Annotation.py --list /data/application/pathongen/data/result/V350071671_20220728163235/78/01.QC/sample.list --type SE --run_mode local --outdir /data/application/pathongen/data/result/V350071671_20220728163235/78/02.Annotation --analysis Taxonomy=PIDB_v1.3.0_rc5,ResistanceGene,VFDB=20210709
################## Upload
/data/Software/Anaconda/anaconda3/bin/python /data/Pipeline/01.PIP-release/PIPv1.1.0rc7/bin/4.Stat/Stat.py --samples /data/application/pathongen/data/result/V350071671_20220728163235/78/sample.list --analysis Merge-QC:Fastp,RemoveHost-Annotation:Taxonomy=PIDB_v1.3.0_rc5,ResistanceGene,VFDB=20210709 --workdir /data/application/pathongen/data/result/V350071671_20220728163235/78 > /data/application/pathongen/data/result/V350071671_20220728163235/78/04.Stat.log 2>&1
################## Waiting for the task to end
wait