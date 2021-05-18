#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(UpSetR)
library(ggplot2)
library(data.table)

FilterFunction <- function(row, svtype) {
  data <- (row["ALT"] == svtype)
}

wanted<-c(1:22,"X","Y")

upset_minimap2<-file.path(args[1],'GM24385.minimap2.upset.tsv')
upsetm_tsv<-fread(upset_minimap2, sep="\t", header=TRUE)
mod<-colnames(upsetm_tsv)
mod[length(mod)]<-"CHROM"
colnames(upsetm_tsv)<-mod
upsetm_tsv$CHROM<-factor(upsetm_tsv$CHROM, levels=wanted)

pdf(file.path(args[1], 'GM24385.minimap2.upset.pdf'), width=15, height=7, onefile=FALSE)
upset(upsetm_tsv, query.legend = "bottom", nintersects = NA, sets = c("TRUTH", "CUTESV", "NPINV", "SNIFFLES", "SVIM", "PBSV"), mb.ratio = c(0.55, 0.45), order.by = "freq", queries = list(list(query = FilterFunction, params = list("DEL"), color = "darkblue", active = F, query.name = "DEL"), list(query = FilterFunction, params = list("INS"), color = "darkred", active = F, query.name = "INS"), list(query = FilterFunction, params = list("DUP"), color = "darkgreen", active = F, query.name = "DUP"), list(query = FilterFunction, params = list("INV"), color = "grey60", active = F, query.name = "INV"),list(query = FilterFunction, params = list("TRA"), color = "yellow", active = F, query.name = "TRA")))
dev.off()

upset_ngmlr<-file.path(args[1],'GM24385.ngmlr.upset.tsv')
upsetn_tsv<-fread(upset_ngmlr, sep="\t", header=TRUE)
mod<-colnames(upsetn_tsv)
mod[length(mod)]<-"CHROM"
colnames(upsetn_tsv)<-mod
upsetn_tsv$CHROM<-factor(upsetn_tsv$CHROM, levels=wanted)

pdf(file.path(args[1], 'GM24385.ngmlr.upset.pdf'), width=15, height=7, onefile=FALSE)
upset(upsetn_tsv, query.legend = "bottom", nintersects = NA, sets = c("TRUTH", "CUTESV", "NPINV", "SNIFFLES", "SVIM"), mb.ratio = c(0.55, 0.45), order.by = "freq", queries = list(list(query = FilterFunction, params = list("DEL"), color = "darkblue", active = F, query.name = "DEL"), list(query = FilterFunction, params = list("INS"), color = "darkred", active = F, query.name = "INS"), list(query = FilterFunction, params = list("DUP"), color = "darkgreen", active = F, query.name = "DUP"), list(query = FilterFunction, params = list("INV"), color = "grey60", active = F, query.name = "INV"),list(query = FilterFunction, params = list("TRA"), color = "yellow", active = F, query.name = "TRA")))
dev.off()

upset_lra<-file.path(args[1],'GM24385.lra.upset.tsv')
upsetn_tsv<-fread(upset_lra, sep="\t", header=TRUE)
mod<-colnames(upsetn_tsv)
mod[length(mod)]<-"CHROM"
colnames(upsetn_tsv)<-mod
upsetn_tsv$CHROM<-factor(upsetn_tsv$CHROM, levels=wanted)

pdf(file.path(args[1], 'GM24385.lra.upset.pdf'), width=15, height=7, onefile=FALSE)
upset(upsetn_tsv, query.legend = "bottom", nintersects = NA, sets = c("TRUTH", "CUTESV", "NPINV", "SNIFFLES", "SVIM"), mb.ratio = c(0.55, 0.45), order.by = "freq", queries = list(list(query = FilterFunction, params = list("DEL"), color = "darkblue", active = F, query.name = "DEL"), list(query = FilterFunction, params = list("INS"), color = "darkred", active = F, query.name = "INS"), list(query = FilterFunction, params = list("INV"), color = "grey60", active = F, query.name = "INV"))) #no DUP/TRA here
dev.off()
