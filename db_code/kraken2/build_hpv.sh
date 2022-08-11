ug_id=1000:1000
docker_name=liumingming1988/biodocker
cwd=`pwd`
db=$cwd/HPV_db
docker_cmd="docker run --rm -v /:/mnt -u $ug_id $docker_name bash -c"
$docker_cmd "kraken2-build --download-taxonomy --db /mnt/$db"
$docker_cmd "kraken2-build --add-to-library /mnt/$cwd/HPV16_18.fa --db /mnt/$db"
$docker_cmd "kraken2-build --build --db /mnt/$db"
