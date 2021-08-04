#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(ggplot2)
library(scales)

fmeasureCurve<- function(f,p) {
  
  return((f*p)/((2*p)-f))
}

plotFMeasures<-function() {
  
  p<-seq(.001,.999,.001)
  f<-seq(.1,.9,.1)
  
  xres<-c()
  yres<-c()
  group<-c()
  
  counter<-0
  gr<-0
  
  for (ftmp in f) {
    gr<-gr+1
    
    for (x in p) {
      
      y<-fmeasureCurve(ftmp, x)
      
      if (0 < y && y <= 1.5) {
        counter<-counter+1
        if (y > 1.03) {
          y <- 1.03
        }
        xres[counter]<-x
        yres[counter]<-y
        group[counter]<-gr
      }
      
    }
    
  }
  
  df<-data.frame(x=xres,y=yres,group=group, stringsAsFactors = FALSE)
  
  return(df)
}


tab<-fread(file.path(args[1],"GM24385.tpfpfn.bylength.tsv"), header=TRUE, sep="\t")

#prbysvtype

result<-list()

aligners<- c("minimap2", "NGMLR", "lra")
tools<-c("cuteSV", "Sniffles", "pbsv", "SVIM")
counter<-0

for (a in aligners) {
  
  for (t in tools) {
    
    #DEL
    
    counter<-counter+1
    stab<-subset(tab, (tool==t & aligner==a & svtype=="DEL"))
    
    tp<-length(which(stab$vartype=="TP"))
    fp<-length(which(stab$vartype=="FP"))
    fn<-length(which(stab$vartype=="FN"))
    
    pdel<-tp/(tp+fp)
    rdel<-tp/(tp+fn)
    
    #INS
    
    stab<-subset(tab, (tool==t & aligner==a & svtype=="INS"))
    
    tp<-length(which(stab$vartype=="TP"))
    fp<-length(which(stab$vartype=="FP"))
    fn<-length(which(stab$vartype=="FN"))

    pins<-tp/(tp+fp)
    rins<-tp/(tp+fn)
    
    r<-data.frame(aligner = a, tool = t, vartype=c("DEL", "INS"), precision=c(pdel,pins), recall=c(rdel,rins))
    result[[counter]]<-r
    
  }
  
}


result<-do.call(rbind,result)
result<-result[complete.cases(result),]

result$aligner<-factor(result$aligner, levels=c("minimap2", "NGMLR", "lra"))
result$vartype<-factor(result$vartype, levels=c("DEL", "INS"))


F1df<-plotFMeasures()
Background<-ggplot() + geom_line(data=F1df,aes(x=x,y=y,group=group),linetype='dashed', color='lightgray') + scale_x_continuous(limits = c(-0.02,1.03), expand=c(0,0), breaks = c(0.0,0.2,0.4,0.6,0.8,1.0)) + scale_y_continuous(expand=c(0,0), breaks=c(0.0,0.2,0.4,0.6,0.8,1.0), limits=c(-0.02,1.03)) + theme_bw() + xlab('Recall') + ylab('Precision') + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + annotate(geom="text", x=0.97, y=0.025, color="black", label='f=0.1', size=3, fontface='italic') + annotate(geom="text", x=0.97, y=0.083, color="black", label='f=0.2', size=3, fontface='italic') + annotate(geom="text", x=0.97, y=0.143, color="black", label='f=0.3', size=3, fontface='italic') + annotate(geom="text", x=0.97, y=0.213, color="black", label='f=0.4', size=3, fontface='italic') + annotate(geom="text", x=0.97, y=0.290, color="black", label='f=0.5', size=3, fontface='italic') + annotate(geom="text", x=0.97, y=0.383, color="black", label='f=0.6', size=3, fontface='italic') + annotate(geom="text", x=0.97, y=0.483, color="black", label='f=0.7', size=3, fontface='italic') + annotate(geom="text", x=0.97, y=0.610, color="black", label='f=0.8', size=3, fontface='italic') + annotate(geom="text", x=0.97, y=0.750, color="black", label='f=0.9', size=3, fontface='italic') + annotate(geom="text", x=0.97, y=0.920, color="black", label='f=1.0', size=3, fontface='italic')
prf1<-Background + geom_point(data=result, aes(x=as.numeric(recall), y=as.numeric(precision), col=tool, shape=aligner), alpha=.8) + facet_wrap(~vartype) +  theme(legend.title=element_blank(), legend.position="bottom", legend.direction = "horizontal",strip.background =element_rect(fill="grey20"), strip.text = element_text(colour = 'white')) +guides(colour = guide_legend(ncol = 7))

ggsave(file.path(args[1],"GM24385.prf1.bysvtype.pdf"), width=12, height=7)


tab$aligner<-factor(tab$aligner, levels=c("minimap2", "NGMLR", "lra"))


p<-ggplot() + geom_histogram(data=tab, aes(size, fill=vartype), binwidth = 0.01)  + theme_bw() + facet_grid(aligner~tool) + labs(x = expression("SV size" [" (log10)"]), y="Count") + theme(legend.title = element_blank(), legend.position = "bottom", legend.direction = "horizontal",strip.background =element_rect(fill="grey20"), strip.text = element_text(colour = 'white')) + scale_x_continuous(trans="log10")


ggsave(file.path(args[1],"GM24385.prf1.bysvlength.pdf"), width=12, height=7)



if (file.exists("Rplots.pdf")) {

    file.remove("Rplots.pdf")
}

