cwd=`pwd`
perl add_pos.pl locate_seq.pos.txt locate_seq.txt|perl get_pos.pl >DJY_predict.bed
sh ../bed/bed_to_interval.sh $cwd/DJY_predict.bed $cwd/DJY_predict.bed.interval_list
