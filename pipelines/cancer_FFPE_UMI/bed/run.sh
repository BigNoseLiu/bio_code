cwd=`pwd`
sh bed_to_interval.sh $cwd/cancer120.bed $cwd/cancer120.bed.interval_list
sh bed_to_cnv_flat.sh $cwd/cancer120.bed $cwd/cancer120.cnvkit.bed
