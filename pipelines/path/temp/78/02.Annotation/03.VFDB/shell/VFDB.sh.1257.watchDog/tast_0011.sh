/data/Software/KMA/KMA-1.3.22/kma -t_db /data/Database/01.PIP-release/db_VFDB/20210709/VFDB -t 8 -ef -i /data/application/pathongen/data/result/V350071671_20220728163235/78/01.QC/02.RemoveHost/N220728902RAAAA/RNA/N220728902RAAAA.fq.gz -apm p -cge -1t1 -o /data/application/pathongen/data/result/V350071671_20220728163235/78/02.Annotation/03.VFDB/N220728902RAAAA/RNA/VFDB
/data/Software/Anaconda/anaconda3/bin/python /data/Pipeline/01.PIP-release/PIPv1.1.0rc7/bin/2.Annotation/3.VFDB/VFDB.filter.py --in /data/application/pathongen/data/result/V350071671_20220728163235/78/02.Annotation/03.VFDB/N220728902RAAAA/RNA/VFDB.res --db 20210709 --identity 95.0 --out /data/application/pathongen/data/result/V350071671_20220728163235/78/02.Annotation/03.VFDB/N220728902RAAAA/RNA/VF.txt
