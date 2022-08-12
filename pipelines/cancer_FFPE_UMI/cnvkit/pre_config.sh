ug_id=1000:1000
docker_name=etal/cnvkit
docker_cmd="docker run --rm -v /:/mnt -u $ug_id $docker_name bash -c"

cwd=`pwd`
bed=$cwd/../bed/cancer120.bed
out_bed=$cwd/cancer120.cnvkit.bed

ref_flat=$cwd/hg19.refFlat.no_chr.txt
ref_fa=$cwd/../../../../databases/gatk_bundle/b37/human_g1k_v37_decoy.fasta
access_bed=$cwd/access-5kb.human_g1k_v37_decoy.bed

#wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz
gunzip -c refFlat.txt.gz|awk '$3!~/_/'|sed 's/chrM/MT/'|sed 's/chr//' >$ref_flat
$docker_cmd "cnvkit.py target /mnt/$bed --annotate /mnt/$ref_flat --split --output /mnt/$out_bed"
#wget https://github.com/etal/cnvkit/blob/master/data/access-5k-mappable.grch37.bed
$docker_cmd "cnvkit.py access /mnt/$ref_fa -s 5000 -o /mnt/$access_bed"
