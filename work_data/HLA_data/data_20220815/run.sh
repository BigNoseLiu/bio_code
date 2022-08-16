#创建目录
cwd=`pwd`
out_dir=$cwd/out_data
if [ ! -d $out_dir ];then
	echo "Success create $out_dir!" 
	date=`date +%Y%M%d%k%M%S`
	real_dir=../../../../../../../../data/analysis_data/out_$date
	mkdir -p $real_dir
	ln -s $real_dir $out_dir
fi

ls /data/raw_data/HLA_data/BFRA228004JX-22803001-T216V1-3s/*gz|perl /home/rain/liumm_data/biodata/git_code/pipelines/HLA_type/AJTK_shell_make.pl >shell_to_run.sh
parallel -j6 --verbose --progress :::: shell_to_run.sh >$out_dir/parallel.log.nohup 2>$out_dir/parallel.log.err
