ug_id=1020:1022
docker_name=liumingming1988/biodocker
biodata_dir=/data/backup/home/liumingming/biodata/
docker_cmd="docker run --rm -v /:/mnt -v $biodata_dir:/biodata --env LD_LIBRARY_PATH=/usr/local/BerkeleyDB/lib/ -u $ug_id $docker_name bash -c "
current_dir=`pwd`
$docker_cmd "python /mnt/$current_dir/test.py >/mnt/$current_dir/out.txt"
