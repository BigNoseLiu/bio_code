#1.如果新建， 去github 新建一个repository

#2. 创建文件夹
git init
git remote add bio_code git@github.com:BigNoseLiu/bio_code.git
git pull bio_code master


#提交改动
#3. 提交改动
git add .
git commit -m 'logs'
git push bio_code master
