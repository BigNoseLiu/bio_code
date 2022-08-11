
#1.如果新建， 去github 新建一个repository

#2.1. 配置秘钥
git config --global user.name "BigNoseLiu"
git config --global user.email "liumingming.2007@163.com"

#2.2. 配置秘钥
ssh-keygen -t rsa -C 'liumingming.2007@163.com'
#将/home/rain/.ssh/id_rsa.pub中的内容拷贝到github的setting/SSH and GPG keys中

#3. 创建文件夹
git init
git remote add bio_code git@github.com:BigNoseLiu/bio_code.git
git pull bio_code master


#4. 提交改动
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

#3. 不同人的修改冲突怎么解决
git stash
git pull
git stash pop	#合并代码
