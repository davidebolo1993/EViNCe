#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(data.table)
library(ComplexUpset)

wanted<-paste0("chr", c(1:22, "X", "Y"))

upset_minimap2<-file.path(args[1],'SI00001.minimap2.upset.tsv')
upsetm_tsv<-fread(upset_minimap2, sep="\t", header=TRUE)
colnames(upsetm_tsv)<-c("ID", "truth", "cuteSV", "npInv", "Sniffles", "SVIM", "pbsv", "ALT", "CHROM")
upsetm_tsv$CHROM<-factor(upsetm_tsv$CHROM, levels=wanted)
upsetm_tsv$ALT<-factor(upsetm_tsv$ALT, levels=c("DEL","INS","DUP","INV","TRA"))

#minimap2

p<-upset(data=upsetm_tsv,intersect=c("truth", "cuteSV", "npInv", "Sniffles", "SVIM", "pbsv"), 
          annotations = list("SV distribution"=(ggplot(mapping=aes(fill=ALT)) + geom_bar(stat='count', position='fill',) + theme_classic() + scale_y_continuous(labels=scales::percent_format()) + ylab('SV distribution') + 
                                                  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(), legend.title=element_blank()))),
          base_annotations=list('Intersection size'=intersection_size(text_colors=c(on_background='black', on_bar='black'), text=list(size=3,angle=45,vjust=0.1,hjust=0.1)) + ylab('Intersection size')),
          #stripes=c('cornsilk1', 'deepskyblue1', 'grey90'),
          width_ratio=0.1,
          keep_empty_groups=TRUE,
          name='Combinations')

ggsave(file.path(args[1], 'SI00001.minimap2.upset.pdf'), width=20, height=15)

#ngmlr

upset_ngmlr<-file.path(args[1],'SI00001.ngmlr.upset.tsv')
upsetm_tsv<-fread(upset_ngmlr, sep="\t", header=TRUE)
colnames(upsetm_tsv)<-c("ID", "truth", "cuteSV", "npInv", "Sniffles", "SVIM", "ALT", "CHROM")
upsetm_tsv$CHROM<-factor(upsetm_tsv$CHROM, levels=wanted)
upsetm_tsv$ALT<-factor(upsetm_tsv$ALT, levels=c("DEL","INS","DUP","INV","TRA"))

p<-upset(data=upsetm_tsv,intersect=c("truth", "cuteSV", "npInv", "Sniffles", "SVIM"), 
          annotations = list("SV distribution"=(ggplot(mapping=aes(fill=ALT)) + geom_bar(stat='count', position='fill',) + theme_classic() + scale_y_continuous(labels=scales::percent_format()) + ylab('SV distribution') + 
                                                  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(), legend.title=element_blank()))),
          base_annotations=list('Intersection size'=intersection_size(text_colors=c(on_background='black', on_bar='black'), text=list(size=3,angle=45,vjust=0.1,hjust=0.1)) + ylab('Intersection size')),
          #stripes=c('cornsilk1', 'deepskyblue1', 'grey90'),
          width_ratio=0.1,
          keep_empty_groups=TRUE,
          name='Combinations')

ggsave(file.path(args[1], 'SI00001.ngmlr.upset.pdf'), width=20, height=15)

#lra

upset_lra<-file.path(args[1],'SI00001.lra.upset.tsv')
upsetm_tsv<-fread(upset_lra, sep="\t", header=TRUE)
colnames(upsetm_tsv)<-c("ID", "truth", "cuteSV", "npInv", "Sniffles", "SVIM", "ALT", "CHROM")
upsetm_tsv$CHROM<-factor(upsetm_tsv$CHROM, levels=wanted)
upsetm_tsv$ALT<-factor(upsetm_tsv$ALT, levels=c("DEL","INS","DUP","INV","TRA"))

p<-upset(data=upsetm_tsv,intersect=c("truth", "cuteSV", "npInv", "Sniffles", "SVIM"), 
          annotations = list("SV distribution"=(ggplot(mapping=aes(fill=ALT)) + geom_bar(stat='count', position='fill',) + theme_classic() + scale_y_continuous(labels=scales::percent_format()) + ylab('SV distribution') + 
                                                  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank(), legend.title=element_blank()))),
          base_annotations=list('Intersection size'=intersection_size(text_colors=c(on_background='black', on_bar='black'), text=list(size=3,angle=45,vjust=0.1,hjust=0.1)) + ylab('Intersection size')),
          #stripes=c('cornsilk1', 'deepskyblue1', 'grey90'),
          width_ratio=0.1,
          keep_empty_groups=TRUE,
          name='Combinations')

ggsave(file.path(args[1], 'SI00001.lra.upset.pdf'), width=20, height=15)


if (file.exists("Rplots.pdf")) {

    file.remove("Rplots.pdf")
}
