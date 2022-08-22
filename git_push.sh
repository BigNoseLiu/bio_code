datetime=`date +%Y%m%d%k%M%S`
git add  --all .
git commit -m "commit at $datetime"
git push bio_code master
