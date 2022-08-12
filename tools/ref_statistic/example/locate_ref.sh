#echo jd A B |perl ../locate_ref.pl -ref ../t1.fa
awk '{print $2"\t"$3}' locate_seq.txt |sed 's/AGATGTGTATAAGAGACAG//g' |perl ../locate_ref.pl -ref ../../../../databases/gatk_bundle/b37/human_g1k_v37_decoy.fasta >locate_seq.result.txt
