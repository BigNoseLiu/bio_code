#创建目录
cwd=`pwd`
out_dir=$cwd/out_data
if [ ! -d $out_dir ];then
	echo 2
	date=`date +%Y%M%d%k%M%S`
	real_dir=../../../../../../../../data/analysis_data/out_$date
	mkdir -p $real_dir
	ln -s $real_dir $out_dir
fi

ls /data/raw_data/DJY/xiaofeng/20220719*gz|perl /home/rain/liumm_data/biodata/git_code/pipelines/DJY_cfDNA_UMI/DJY_shell_make.pl  >shell_to_run.sh
nohup parallel -j6 --verbose --progress :::: shell_to_run.sh >$out_dir/nohup 2>$out_dir/err &
