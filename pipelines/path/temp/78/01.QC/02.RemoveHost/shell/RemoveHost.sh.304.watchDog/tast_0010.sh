/data/Software/SNAP/SNAP-v1.0.0/snap-aligner single /data/Database/01.PIP-release/db_Host/SNAP/20220620 -t 14 -compressedFastq /data/application/pathongen/data/result/V350071671_20220728163235/78/01.QC/01.Fastp/N220728903HAAAA/DNA/N220728903HAAAA.fq.gz -o -bam - |/data/Software/Anaconda/anaconda3/envs/Samtools/bin/samtools sort -@ 14 -o /data/application/pathongen/data/result/V350071671_20220728163235/78/01.QC/02.RemoveHost/N220728903HAAAA/DNA/N220728903HAAAA.bam -
/data/Software/bamdst/bamdst-v1.0.9/bamdst -q 1 -p /data/Database/01.PIP-release/db_Host/SNAP/20220620/InterTarget.bed -o /data/application/pathongen/data/result/V350071671_20220728163235/78/01.QC/02.RemoveHost/N220728903HAAAA/DNA /data/application/pathongen/data/result/V350071671_20220728163235/78/01.QC/02.RemoveHost/N220728903HAAAA/DNA/N220728903HAAAA.bam &
/data/Software/Anaconda/anaconda3/envs/Samtools/bin/samtools fastq -f 4 -@ 13 /data/application/pathongen/data/result/V350071671_20220728163235/78/01.QC/02.RemoveHost/N220728903HAAAA/DNA/N220728903HAAAA.bam -0 /data/application/pathongen/data/result/V350071671_20220728163235/78/01.QC/02.RemoveHost/N220728903HAAAA/DNA/N220728903HAAAA.fq.gz &
wait
