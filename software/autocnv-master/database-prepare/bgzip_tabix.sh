type=$1
file=$2
bgzip -cf $file >$file.gz
tabix -fp $type $file.gz
rm -rf $file
