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
parallel -j6 --verbose --progress :::: shell_to_run.sh >$out_dir/parallel.log.nohup 2>$out_dir/parallel.log.err

#stat vcf
ls out_data/*/*hg19_multianno.vcf|perl ../../../pipelines/cancer_FFPE_UMI/merge_vcf.v1.02.pl |perl ../../../pipelines/DJY_cfDNA_UMI/DJY_convert_vcf.v1.01.pl >all_vcf.merge.xls

#stat qc
ls out_data/*/*.08.*metrics out_data/*/*.05.*metrics |perl ../../../pipelines/cancer_FFPE_UMI/stat_picardQC.v1.01.pl recal.qc_stat.xls
ls out_data/*/*.07.*metrics out_data/*/*.05.*metrics |perl ../../../pipelines/cancer_FFPE_UMI/stat_picardQC.v1.01.pl sorted.qc_stat.xls
