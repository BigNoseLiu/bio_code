ug_id=1020:1022
docker_name=liumingming1988/biodocker
biodata_dir=/data/backup/home/liumingming/biodata/
docker_cmd="docker run --rm -v /:/mnt -v $biodata_dir:/biodata --env LD_LIBRARY_PATH=/usr/local/BerkeleyDB/lib/ -u $ug_id $docker_name bash -c "
current_dir=`pwd`


mkdir -p out_target out_stat

ls raw_panels|sed 's/.gene_list//'|while read panel_name
do
	cat raw_panels/$panel_name\.gene_list|perl ../bin/get_target_gff.pl /data/backup/home/liumingming/biodata/databases/gene_info/GCF_000001405.25_GRCh37.p13_genomic.gff.gz out_target/$panel_name


	ls out_target/*/*/*.CDS out_target/*/*/*.exon |while read L
	do
		cat $L|sort -k 1,1 -k 2,2n -k 3,3n >$L\.sorted.bed
		$docker_cmd "bedtools merge -i /mnt/$current_dir/$L\.sorted.bed -c 4,5 -o distinct,distinct -delim \"|\" >/mnt/$current_dir/$L\.sorted.merge.bed"
	done


	for type in 'CDS' 'exon'
	do
		cat out_target/$panel_name/$type/*sorted.merge.bed |sort -k 1,1 -k 2,2n -k 3,3n >out_stat/$panel_name\.$type\.merge.bed
		ls target_bed|while read panel
		do
			$docker_cmd "bedtools coverage -a  /mnt/$current_dir/out_stat/$panel_name\.$type\.merge.bed  -b /mnt/$current_dir/target_bed/$panel >/mnt/$current_dir/out_stat/$panel_name\.$type\_cover_$panel\.cover.xls"
		done
	done
done

ls out_stat/*cover.xls|perl ../bin/stat_coverage.pl >out_stat/gene_cover.stat.xls 2>out_stat/panel_cover.stat.xls

