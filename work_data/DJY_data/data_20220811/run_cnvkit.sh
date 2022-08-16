#创建目录
cwd=`pwd`
out_dir=$cwd/cnv_out
if [ ! -d $out_dir ];then
	echo "Success create $out_dir!" 
	date=`date +%Y%M%d%k%M%S`
	real_dir=../../../../../../../../data/analysis_data/out_$date
	mkdir -p $real_dir
	ln -s $real_dir $out_dir
fi
sh /home/rain/liumm_data/biodata/git_code/pipelines/DJY_cfDNA_UMI/cnvkit/cnvkit_noRef.sh 1000:1000 cnv_cell $out_dir  $cwd/cnvkit_in/dis_cell/ $cwd/cnvkit_in/normal_cell
sh /home/rain/liumm_data/biodata/git_code/pipelines/DJY_cfDNA_UMI/cnvkit/cnvkit_noRef.sh 1000:1000 cnv_cfDNA $out_dir $cwd/cnvkit_in/dis_cfDNA/ $cwd/cnvkit_in/normal_cfDNA
