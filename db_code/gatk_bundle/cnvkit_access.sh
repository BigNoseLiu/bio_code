
ug_id=1000:1000
docker_name=etal/cnvkit
docker_cmd="docker run --rm -v /:/mnt -u $ug_id $docker_name bash -c"
cwd=`pwd`
$docker_cmd "cnvkit.py access /mnt/$cwd/human_g1k_v37_decoy.fasta -s 5000 -o /mnt/$cwd/access-5kb.human_g1k_v37_decoy.bed"
