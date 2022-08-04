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




#FAQ
#1. 遇到超过>100M的文件，移除相关的文件，修改git的历史记录，移除相应的commit结点。然后才能push
git filter-branch -f --index-filter 'git rm --cached --ignore-unmatch software/annovar/annovar.latest.tar.gz'

#2. github如果访问不到，
修改/etc/hosts文件:  sudo vi /etc/hosts
添加如下行
140.82.113.3      github.com
