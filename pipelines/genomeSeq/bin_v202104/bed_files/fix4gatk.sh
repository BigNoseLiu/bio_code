bed_file=$1
cdir=`pwd`;
cat $bed_file\.bed|perl fix_bed_format_for_gatk.pl >$bed_file.fix4gatk.bed
java -jar ../../../software/Picard/20180930/picard.jar BedToIntervalList I=$bed_file.fix4gatk.bed O=$bed_file.fix4gatk.bed.interval_list SD=../../../database/genomes/b37_decoy/human_g1k_v37_decoy.dict
docker run -w /out_dir -u 1002:1002 --rm -v /home/liumingming/bio_pipeline/dockers/biodocker/biodata:/biodata -v :/out_dir li/biojupyter bash -c "  picard
