#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(RColorBrewer)
library(ggplot2)
library(viridis)
library(dplyr)
library(gridExtra)

#plot by chromosome

wanted<-c(1:22,"X","Y")

#minimap2 plot
minimap2plots<-list()
tab<-fread(file.path(args[1], "GM24385.minimap2.regions.bed.gz"), sep="\t")
#remove unwanted chromosomes and exclude biggest outlyers
forylab<-c("1","5","9","13","17","21")
forxlab<-c("21", "22", "X", "Y")
for (i in 1:length(wanted)) {
  #proceed by chromosome
  subtab<-data.frame(subset(tab, (V1 == wanted[i])), stringsAsFactors = FALSE)
  m<-mean(subtab$V4)
  d<-sd(subtab$V4)
  subtab$V5<-(subtab$V4-m)/d
  zlimit<-1.5
  subsubtab<-data.frame(subset(subtab, (V5 <= zlimit)), stringsAsFactors = FALSE)
  excluded<-data.frame(subset(subtab, (V5 > zlimit)), stringsAsFactors = FALSE)
  meand<-mean(subsubtab$V4)
  downsampled<-subsubtab %>% sample_frac(.5) #plot a subsample o
  ylab_<-""
  xlab_<-""
  if (wanted[i] %in% forylab) {
    ylab_<-"Depth"
  }
  if (wanted[i] %in% forxlab) {
    xlab_<-"Genomic coordinate"
  }
  minimap2plots[[wanted[i]]]<-ggplot(downsampled, aes(x=V2, y=V4, color=V4)) + geom_point(size=.2, show.legend = FALSE) + scale_color_viridis("Depth", option = "plasma") + theme_bw() + xlab(xlab_) + ylab(ylab_) + scale_y_continuous(limits=c(0, max(downsampled$V4)+20)) + geom_hline(yintercept=meand, linetype="dashed", color = "black") + ggtitle(paste0("Chromosome ", wanted[i])) + theme(plot.title = element_text(hjust = 0.5))
  fwrite(excluded, file=file.path(args[2], "GM24385.minimap2.exclude.tsv"), sep="\t",quote=FALSE, col.names=FALSE, row.names = FALSE, append=TRUE)
}
  
n <- length(minimap2plots)
ncols <- floor(sqrt(n))
global<-do.call("grid.arrange", c(minimap2plots, ncol=ncols))

ggsave(file.path(args[1], "GM24385.minimap2.coverage_per_chromosome.png"), global, width=20, height = 10)

#ngmlr plot
ngmlrplots<-list()
tab<-fread(file.path(args[1], "GM24385.ngmlr.regions.bed.gz"), sep="\t")
#remove unwanted chromosomes and exclude biggest outlyers
for (i in 1:length(wanted)) {
  #proceed by chromosome
  subtab<-data.frame(subset(tab, (V1 == wanted[i])), stringsAsFactors = FALSE)
  m<-mean(subtab$V4)
  d<-sd(subtab$V4)
  subtab$V5<-(subtab$V4-m)/d
  zlimit<-1.5
  subsubtab<-data.frame(subset(subtab, (V5 <= zlimit)), stringsAsFactors = FALSE)
  excluded<-data.frame(subset(subtab, (V5 > zlimit)), stringsAsFactors = FALSE)
  meand<-mean(subsubtab$V4)
  downsampled<-subsubtab %>% sample_frac(.5)
  ylab_<-""
  xlab_<-""
  if (wanted[i] %in% forylab) {
    ylab_<-"Depth"
  }
  if (wanted[i] %in% forxlab) {
    xlab_<-"Genomic coordinate"
  }
  ngmlrplots[[wanted[i]]]<-ggplot(downsampled, aes(x=V2, y=V4, color=V4)) + geom_point(size=.2, show.legend = FALSE) + scale_color_viridis("Depth", option = "plasma") + theme_bw() + xlab(xlab_) + ylab(ylab_) + scale_y_continuous(limits=c(0, max(downsampled$V4)+20)) + geom_hline(yintercept=meand, linetype="dashed", color = "black") + ggtitle(paste0("Chromosome ", wanted[i])) + theme(plot.title = element_text(hjust = 0.5))
  fwrite(excluded, file=file.path(args[2], "GM24385.ngmlr.exclude.tsv"), sep="\t",quote=FALSE, col.names=FALSE, row.names = FALSE, append=TRUE)
}

n <- length(ngmlrplots)
ncols <- floor(sqrt(n))
global<-do.call("grid.arrange", c(ngmlrplots, ncol=ncols))
ggsave(file.path(args[1], "GM24385.ngmlr.coverage_per_chromosome.png"), global, width=20, height = 10)

#pbmm2 plot
pbmm2plots<-list()
tab<-fread(file.path(args[1], "GM24385.pbmm2.regions.bed.gz"), sep="\t")
#remove unwanted chromosomes and exclude biggest outlyers
for (i in 1:length(wanted)) {
  #proceed by chromosome
  subtab<-data.frame(subset(tab, (V1 == wanted[i])), stringsAsFactors = FALSE)
  m<-mean(subtab$V4)
  d<-sd(subtab$V4)
  subtab$V5<-(subtab$V4-m)/d
  zlimit<-1.5
  subsubtab<-data.frame(subset(subtab, (V5 <= zlimit)), stringsAsFactors = FALSE)
  excluded<-data.frame(subset(subtab, (V5 > zlimit)), stringsAsFactors = FALSE)
  meand<-mean(subsubtab$V4)
  downsampled<-subsubtab %>% sample_frac(.5)
  ylab_<-""
  xlab_<-""
  if (wanted[i] %in% forylab) {
    ylab_<-"Depth"
  }
  if (wanted[i] %in% forxlab) {
    xlab_<-"Genomic coordinate"
  }
  pbmm2plots[[wanted[i]]]<-ggplot(downsampled, aes(x=V2, y=V4, color=V4)) + geom_point(size=.2, show.legend = FALSE) + scale_color_viridis("Depth", option = "plasma") + theme_bw() + xlab(xlab_) + ylab(ylab_) + scale_y_continuous(limits=c(0, max(downsampled$V4)+20)) + geom_hline(yintercept=meand, linetype="dashed", color = "black") + ggtitle(paste0("Chromosome ", wanted[i])) + theme(plot.title = element_text(hjust = 0.5))
  fwrite(excluded, file=file.path(args[2], "GM24385.pbmm2.exclude.tsv"), sep="\t",quote=FALSE, col.names=FALSE, row.names = FALSE, append=TRUE)
}

n <- length(pbmm2plots)
ncols <- floor(sqrt(n))
global<-do.call("grid.arrange", c(pbmm2plots, ncol=ncols))
ggsave(file.path(args[1], "GM24385.pbmm2.coverage_per_chromosome.png"), global, width=20, height = 10)

#lra plot
lraplots<-list()
tab<-fread(file.path(args[1], "GM24385.lra.regions.bed.gz"), sep="\t")
#remove unwanted chromosomes and exclude biggest outlyers
for (i in 1:length(wanted)) {
  #proceed by chromosome
  subtab<-data.frame(subset(tab, (V1 == wanted[i])), stringsAsFactors = FALSE)
  m<-mean(subtab$V4)
  d<-sd(subtab$V4)
  subtab$V5<-(subtab$V4-m)/d
  zlimit<-1.5
  subsubtab<-data.frame(subset(subtab, (V5 <= zlimit)), stringsAsFactors = FALSE)
  excluded<-data.frame(subset(subtab, (V5 > zlimit)), stringsAsFactors = FALSE)
  meand<-mean(subsubtab$V4)
  downsampled<-subsubtab %>% sample_frac(.5)
  ylab_<-""
  xlab_<-""
  if (wanted[i] %in% forylab) {
    ylab_<-"Depth"
  }
  if (wanted[i] %in% forxlab) {
    xlab_<-"Genomic coordinate"
  }
  lraplots[[wanted[i]]]<-ggplot(downsampled, aes(x=V2, y=V4, color=V4)) + geom_point(size=.2, show.legend = FALSE) + scale_color_viridis("Depth", option = "plasma") + theme_bw() + xlab(xlab_) + ylab(ylab_) + scale_y_continuous(limits=c(0, max(downsampled$V4)+20)) + geom_hline(yintercept=meand, linetype="dashed", color = "black") + ggtitle(paste0("Chromosome ", wanted[i])) + theme(plot.title = element_text(hjust = 0.5))
  fwrite(excluded, file=file.path(args[2], "GM24385.lra.exclude.tsv"), sep="\t",quote=FALSE, col.names=FALSE, row.names = FALSE, append=TRUE)
}

n <- length(lraplots)
ncols <- floor(sqrt(n))
global<-do.call("grid.arrange", c(lraplots, ncol=ncols))
ggsave(file.path(args[1], "GM24385.lra.coverage_per_chromosome.png"), global, width=20, height = 10)


#minimap2,ngmlr,pbmm2 and lra plot, fraction of bases per depth treshold
files<-c(file.path(args[1], "GM24385.minimap2.mosdepth.global.dist.txt"), file.path(args[1], "GM24385.ngmlr.mosdepth.global.dist.txt"), file.path(args[1],"GM24385.pbmm2.mosdepth.global.dist.txt"),file.path(args[1],"GM24385.lra.mosdepth.global.dist.txt"))

png(file.path(args[1], "GM24385.fraction_of_bases_per_depth.png"), height=1000, width=1000)
plot(c(1:200), c(1:200), type='n', xlab="Depth", ylab="Fraction of bases \u2265 depth", ylim=c(0,1.0), font.main = 1, main="Proportion of genome at coverage")
abline(v = 20, col = "gray60")
abline(v = 50, col = "gray60")
abline(v = 80, col = "gray60")
abline(v = 100, col = "gray60")
abline(v = 200, col = "gray60")
abline(v = 300, col = "gray60")
abline(v = 400, col = "gray60")
abline(h = 0.50, col = "gray60")
abline(h = 0.90, col = "gray60")
axis(1, at=c(20,50,80), labels=c(20,50,80))
axis(2, at=c(0.50), labels=c(0.50))

labels_<-c("minimap2", "ngmlr", "pbmm2", "lra")
colors_<-c("darkred", "darkblue", "darkgreen", "black")
tab<-fread(files[1],sep="\t")
subsubtab<-subset(tab, (V1 == "total" & V2 < 200))
points(subsubtab$V2, subsubtab$V3, type='l', lwd=3, col=colors_[1])
tab<-fread(files[2],sep="\t")
subsubtab<-subset(tab, (V1 == "total" & V2 < 200))
points(subsubtab$V2, subsubtab$V3, type='l', lwd=3, col=colors_[2])
tab<-fread(files[3],sep="\t")
subsubtab<-subset(tab, (V1 == "total" & V2 < 200))
points(subsubtab$V2, subsubtab$V3, type='l', lwd=3, col=colors_[3])
tab<-fread(files[4],sep="\t")
subsubtab<-subset(tab, (V1 == "total" & V2 < 200))
points(subsubtab$V2, subsubtab$V3, type='l', lwd=3, col=colors_[4])

legend("topright", legend=labels_, col=colors_, lty=1, lwd=4)
dev.off()

if (file.exists("Rplots.pdf")) {

    file.remove("Rplots.pdf")
}



