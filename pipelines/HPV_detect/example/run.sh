cwd=`pwd`
mkdir -p out_dir
ls /data/raw_data/genolab/202207131614*/*/*gz|grep -v NA|grep -v phix|perl /home/rain/liumm_data/biodata/pipelines/HPV_detect/shell_make.pl $cwd/out_dir >all_shell.sh
ls /data/raw_data/genolab/202207011442_220101009_2P2206040*/*/*gz|grep -v NA|grep -v phix|perl /home/rain/liumm_data/biodata/pipelines/HPV_detect/shell_make.pl $cwd/out_dir >>all_shell.sh
ls /data/raw_data/genolab/*/*/*gz|grep HY|perl /home/rain/liumm_data/biodata/pipelines/HPV_detect/shell_make.pl $cwd/out_dir >>all_shell.sh
ls out_dir/POOLINHG*/*report|perl /home/rain/liumm_data/biodata/pipelines/HPV_detect/stat_hpv.pl >hpv_220718.xls
