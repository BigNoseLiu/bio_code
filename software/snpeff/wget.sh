#http://snpeff.sourceforge.net/download.html#databases
#1. download
#2. config the  snpEff/snpEff.config file
#3. download database from https://sourceforge.net/projects/snpeff/files/databases/v4_3/
unzip snpEff_v4_3_GRCh37.87.zip

#或者直接下载
java -jar snpEff.jar download GRCh37.p13.RefSeq

#4. 移动data到database目录
#5. edit snpEff/snpEff.config,添加如下行 设定数据目录
data.dir = ../../../../databases/snpeff_data/

