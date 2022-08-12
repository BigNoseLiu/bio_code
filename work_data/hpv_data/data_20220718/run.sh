#创建目录
cwd=`pwd`
out_dir=$cwd/out_data
if [ ! -d $out_dir ];then
	date=`date +%Y%M%d%k%M%S`
	real_dir=../../../../../../../../data/analysis_data/out_$date
	mkdir -p $real_dir
	ln -s $real_dir $out_dir
	echo "Success create $out_dir!" 
fi

rm -rf shell_to_run.sh
for tag in 'HY' 'ZMH2Y' 'ZMHBY'
do
	find /data/raw_data/genolab/ -type f|grep $tag|perl /home/rain/liumm_data/biodata/git_code/pipelines/HPV_detect/shell_make.pl >>shell_to_run.sh
done
parallel -j2 --verbose --progress :::: shell_to_run.sh >$out_dir/parallel.log.nohup 2>$out_dir/parallel.log.err

#stat all
ls out_data/*/*report|perl ../../../pipelines/HPV_detect/stat_hpv.pl >hpv_stat.xls
