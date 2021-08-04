#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(ggplot2)

files<-c(file.path(args[1], "GM24385.minimap2.mosdepth.global.dist.txt"), file.path(args[1], "GM24385.ngmlr.mosdepth.global.dist.txt"), file.path(args[1],"GM24385.pbmm2.mosdepth.global.dist.txt"),file.path(args[1],"GM24385.lra.mosdepth.global.dist.txt"))

tab<-fread(files[1], sep="\t")
tab1<-data.frame(subset(tab, (V1 == "total" & V2 < 200)))
colnames(tab1)<-c("chr", "cov", "perc")
tab1$aligner<-"minimap2"

tab<-fread(files[2], sep="\t")
tab2<-data.frame(subset(tab, (V1 == "total" & V2 < 200)))
colnames(tab2)<-c("chr", "cov", "perc")
tab2$aligner<-"NGMLR"

tab<-fread(files[3], sep="\t")
tab3<-data.frame(subset(tab, (V1 == "total" & V2 < 200)))
colnames(tab3)<-c("chr", "cov", "perc")
tab3$aligner<-"pbmm2"

tab<-fread(files[4], sep="\t")
tab4<-data.frame(subset(tab, (V1 == "total" & V2 < 200)))
colnames(tab4)<-c("chr", "cov", "perc")
tab4$aligner<-"lra"

taball<-rbind(tab1,tab2,tab3,tab4)
taball$aligner<-factor(taball$aligner, levels=c("minimap2", "NGMLR", "pbmm2", "lra"))

p<-ggplot(taball, aes(x=cov, y=perc,col=aligner)) + geom_line() + labs(x="Depth", y=expression("Fraction of bases ">=" depth")) + theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom", legend.direction = "horizontal") + scale_x_continuous(limits=c(0,200))

ggsave(file.path(args[1], "GM24385.fraction_of_bases_per_depth.pdf"), device="pdf", height=10, width=7)

if (file.exists("Rplots.pdf")) {

    file.remove("Rplots.pdf")
}




