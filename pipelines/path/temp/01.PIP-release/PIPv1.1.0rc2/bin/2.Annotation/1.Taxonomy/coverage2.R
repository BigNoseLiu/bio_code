args = commandArgs(T)
library(ggplot2)
dat = read.table(args[1], sep='\t')
colnames(dat) = c('latin', 'bin', 'depth')
for (nc in unique(dat$latin))
{
    dat2 = dat[which(dat$latin==nc),]
    sps = gsub("/", "_", nc)
    out = paste(args[2], paste(gsub(" ","_",sps), "coverage.png", sep="."), sep="/")
    out2 = gsub("\\(", "\\\\(", out)
    out3 = gsub("\\)", "\\\\)", out2)
    png(out3, width=900, height=500)
    #bitmap(paste(args[2], paste(gsub("/", "_", gsub(" ","_",nc)), "coverage.png", sep="."), sep="/"), width=900, height=500,"png16m")
    print(ggplot(data=dat2, aes(x=bin, y=depth))+geom_bar(stat="identity")+xlab("Position")+ylab("Depth")+ggtitle(paste(nc, " coverage=", round(100*nrow(dat2[which(dat2$depth>0),])/nrow(dat2), 2), "%", sep=""))+theme_bw(base_size=20, base_family="DejaVu Sans")+guides(fill=FALSE) + theme(plot.title = element_text(hjust = 0.5)))
    dev.off()
}
