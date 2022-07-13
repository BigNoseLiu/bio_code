datetime=`date -d@1234567890 +"%F %T"`
git add  --all .
git commit -m "commit $datetime"
git push bio_code master
