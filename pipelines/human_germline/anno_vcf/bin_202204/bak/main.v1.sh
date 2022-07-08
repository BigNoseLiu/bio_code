#docker settings
ug_id=1020:1022
docker_name=liumingming1988/biodocker
biodata_dir=/data/backup/home/liumingming/biodata/
docker_cmd="docker run --rm -v /:/mnt -v $biodata_dir:/biodata --env LD_LIBRARY_PATH=/usr/local/BerkeleyDB/lib/ -u $ug_id $docker_name bash -c "
current_dir=`pwd`


#input
sample_id=$1
out_dir=$2/$sample_id
vcf_in=$3
mkdir -p $out_dir



#sofewares
bin_snpeff="java -jar /biodata/software/snpeff/snpEff5_0/snpEff.jar -v GRCh37.p13.RefSeq";
#bin_select_snpeff	=	"$docker_cmd perl $bin_dir/snpeff/select_standard_hgvs_from_snpeff_Eff.v1.01.pl";
#bin_fix_snpeff	=	"$docker_cmd perl $bin_dir/snpeff/fix_snpeff_Eff.v1.pl";


#out files



#commnads
