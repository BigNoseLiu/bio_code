#在线安装包
#进入R界面
#R-3.5.1/bin/R
#install.packages(/share/liumingming/XArchive/bio_pipeline/software/R/lib_gz/cghFLasso_0.2-1.tar.gz, repos = NULL, type="source")
#install.packages(/share/liumingming/XArchive/bio_pipeline/software/R/lib_gz/reshape_0.8.8.tar.gz, repos = NULL, type="source")

#运行脚本
../../software/R/R-3.5.1/bin/Rscript  pedigree.R kinship.csv outplot_pedigree --no-save
