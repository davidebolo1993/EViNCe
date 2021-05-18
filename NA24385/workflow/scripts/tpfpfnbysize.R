#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(ggplot2)

tab<-fread(file.path(args[1],"GM24385.tpfpfn.bylength.tsv"), header=TRUE, sep="\t")
colors<-c("TP" = "#440154FF", "FP"="#21908CFF", "FN" =  "#FDE725FF", "DEL" ="#0D0887FF", "INS" = "#ED7953FF")
p<-ggplot() + geom_histogram(data=tab, aes(size, fill=vartype), binwidth = 0.01)  + theme_bw() + facet_grid(aligner~tool) + labs(x = expression("SV size" [" (log10)"]), y="Count") + theme(legend.title = element_blank(), legend.position = "bottom", legend.direction = "horizontal",strip.background =element_rect(fill="grey20"), strip.text = element_text(colour = 'white')) + scale_fill_manual(values=colors,breaks=c("DEL","INS", "FP", "TP", "FN")) + scale_x_continuous(trans="log10")
ggsave(file.path(args[1],"GM24385.tpfpfn.bylength.pdf"), width=12, height=7) 
p<-ggplot() + geom_histogram(data=tab, aes(size, fill=svtype), binwidth = 0.01)  + theme_bw() + facet_grid(aligner~tool) + labs(x = expression("SV size" [" (log10)"]), y="Count") + theme(legend.title = element_blank(), legend.position = "bottom", legend.direction = "horizontal",strip.background =element_rect(fill="grey20"), strip.text = element_text(colour = 'white')) + scale_fill_manual(values=colors,breaks=c("DEL","INS", "FP", "TP", "FN")) + scale_x_continuous(trans="log10")
ggsave(file.path(args[1],"GM24385.svtype.bylength.pdf"), width=12, height=7) 


if (file.exists("Rplots.pdf")) {

    file.remove("Rplots.pdf")
}

