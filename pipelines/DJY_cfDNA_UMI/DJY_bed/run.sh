cwd=`pwd`
perl add_pos.pl locate_seq.result.txt locate_seq.txt|perl get_pos.v1.02.pl|sort -k 1,1 -k 2,3n >DJY_predict.bed
sh $cwd/../../cancer_FFPE_UMI/bed/bed_to_interval.sh $cwd/DJY_predict.bed $cwd/DJY_predict.bed.interval_list
