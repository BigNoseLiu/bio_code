docker run --rm -v /:/mnt -v /data/backup/home/liumingming/biodata:/biodata --env LD_LIBRARY_PATH=/usr/local/BerkeleyDB/lib/ -u 1020:1022 liumingming1988/biodocker  bash -c "cat /mnt//data/backup/home/liumingming/biodata/git_code/pipelines/human_germline/anno_vcf/example/out_data/merge_result/merge_all.final.leftalign.annoGene.annoVar.hg19_merge.filter.txt|perl /biodata/git_code/pipelines/human_germline/anno_vcf/bin_202204/mongodb/anno2mongodb.v1.01.pl Genetic.Tests merge_all "
