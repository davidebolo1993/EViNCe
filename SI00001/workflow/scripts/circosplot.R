#!/usr/bin/env Rscript

library(grid)
library(circlize)
library(data.table)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ComplexHeatmap)

#legend 

lgd_points = Legend(at = c("INS", "DEL", "DUP", "INV"), type = "points", 
    legend_gp = gpar(col = c("#39558CFF", "#440154FF", "#FDE725FF", "#74D055FF")), title_position = "topleft", 
    title = "", nrow = 1, background="#FFFFFF")

lgd_lines = Legend(at = c("TRA"), type = "lines", 
    legend_gp = gpar(col = "black" , lwd = 2), title_position = "topleft", 
    title = "", nrow = 1,background="#FFFFFF")

lgd_list_horizontal = packLegend(lgd_points, lgd_lines, direction = "horizontal")
circle_size = unit(1, "snpc")

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 5) {

  tool<-args[4]
  outdir<-file.path(args[5])
  pdf(file.path(outdir, paste0("SI00001.", tool, ".svs.pdf")), height=6, width=15)
  layout(matrix(1:3, 1,3))

  for (dir in c(args[1], args[2], args[3])) {

    bedpe<-file.path(dir, "SI00001.bnd.bedpe") #this exists, but can be empty
    info<-file.info(bedpe)

    if (info$size != 0) { #this is not always true

      BNDs<-fread(bedpe, sep="\t", header=FALSE)
      colnames(BNDs)<-c("chr", "pos1", "end1", "chr2", "pos2", "end2", "type", "num", "s1", "s2")

    }
          
    bed<-file.path(dir, "SI00001.nobnd.bed")
    info2<-file.info(bed)

    if (info2$size == 0) {
      #empty ploy
      circos.par("start.degree" = 90, circle.margin=0.0001)
      circos.initializeWithIdeogram(species = "hg38")
      next
    }

    NOBNDs<-fread(bed, sep="\t", header=FALSE)
    colnames(NOBNDs)<-c("chr", "pos", "end", "type")
    NOBNDRanges<-makeGRangesFromDataFrame(NOBNDs, keep.extra.columns = TRUE,start.field = "pos", end.field = "end")
    chrSizes<-seqlengths(Hsapiens)[c(1:24)]
    windows<-tileGenome(chrSizes, tilewidth=1e7, cut.last.tile.in.chrom=T)
    svtypes<-c("INS", "DEL", "DUP", "INV")
    finvalues<-list()
          
    for (l in 1:length(windows)) {
            
      w<-windows[l]
      chr<-as.character(w@seqnames@values)
      start<-as.numeric(w@ranges@start)
      end<-as.numeric(data.frame(w@ranges)$end)
      bps<-as.numeric(w@ranges@width)
      res<-findOverlaps(w, NOBNDRanges)
      idx<-res@to
      subw<-NOBNDs[idx,]
            
      #iterate over types
            
      values<-list()
            
      for (i in 1:length(svtypes)) {
              
        idxtype<-which(subw$type == svtypes[i])
              
        if (length(idxtype) == 0) {
                
          ratio<-0 #no alt
                
        } else {
                
          subsubw<-subw[idxtype,]      
          diff<-subsubw$end-subsubw$pos
          cums<-cumsum(diff)[length(diff)]
          ratio<-cums/bps
                
        }
              
        values[[svtypes[i]]]<-ratio
              
      }
            
      block<-do.call(c, values)
      vector<-data.frame('chr' = chr, 'start' = start, 'end' = end, 'ins' = as.numeric(block[1]), 'del' = as.numeric(block[2]), 'dup' = as.numeric(block[3]), 'inv' = as.numeric(block[4]))
      finvalues[[l]]<-vector
            
    }
          
    heattable<-do.call(rbind,finvalues)
    heatdataframe<-split(data.table(heattable), by="chr")
    reslist<-list()

    for (el in names(heatdataframe)) {
          
      if (min(heatdataframe[[el]]$ins) == max(heatdataframe[[el]]$ins)) {
              
        inscolors<-rep("#FFFFFF",nrow(heatdataframe[[el]]))
              
      } else {
            
        ii<-cut(heatdataframe[[el]]$ins, breaks = seq(min(heatdataframe[[el]]$ins), max(heatdataframe[[el]]$ins), len = 100), include.lowest = TRUE)
        inscolors <- colorRampPalette(c("white", "darkblue"))(99)[ii]
              
      }

      if (min(heatdataframe[[el]]$del) == max(heatdataframe[[el]]$del)) {
              
        delcolors<-rep("#FFFFFF",nrow(heatdataframe[[el]]))
              
      } else {
              
        ii<-cut(heatdataframe[[el]]$del, breaks = seq(min(heatdataframe[[el]]$del), max(heatdataframe[[el]]$del), len = 100), include.lowest = TRUE)
        delcolors <- colorRampPalette(c("white", "mediumorchid4"))(99)[ii]
              
      }
            
      if (min(heatdataframe[[el]]$dup) == max(heatdataframe[[el]]$dup)) {
              
        dupcolors<-rep("#FFFFFF",nrow(heatdataframe[[el]]))
              
      } else {
              
        ii<-cut(heatdataframe[[el]]$dup, breaks = seq(min(heatdataframe[[el]]$dup), max(heatdataframe[[el]]$dup), len = 100), include.lowest = TRUE)
        dupcolors <- colorRampPalette(c("white", "gold1"))(99)[ii]
              
      }
            
      if (min(heatdataframe[[el]]$inv) == max(heatdataframe[[el]]$inv)) {
              
        invcolors<-rep("#FFFFFF",nrow(heatdataframe[[el]]))
              
      } else {
              
        ii<-cut(heatdataframe[[el]]$inv, breaks = seq(min(heatdataframe[[el]]$inv), max(heatdataframe[[el]]$inv), len = 100), include.lowest = TRUE)
        invcolors <- colorRampPalette(c("white", "darkgreen"))(99)[ii]
              
      }
            
      df<-cbind(inscolors,delcolors,dupcolors,invcolors)
      reslist[[el]]<-df
    
    }
         
    mat_col<-matrix(do.call(rbind,reslist),ncol = 4)
          
    circos.par("start.degree" = 90, circle.margin=0.0001)
    circos.initializeWithIdeogram(species = "hg38")
    circos.genomicHeatmap(heattable, col = mat_col,line_col = as.numeric(factor(heatdataframe[[1]])),connection_height = NULL)
          
    NONBNSs_colors<-list()
    NONBNSs_colors[["INS"]]<-"#39558CFF"
    NONBNSs_colors[["DEL"]]<-"#440154FF"
    NONBNSs_colors[["DUP"]]<-"#FDE725FF"
    NONBNSs_colors[["INV"]]<-"#74D055FF"
          
    #PLOT NON-BNDs first
          
    bed_list<-split(NOBNDs,by="type")
          
    for (el in svtypes) {
            
      if (el %in% names(bed_list)) {
              
        circos.genomicDensity(bed_list[[el]], col = NONBNSs_colors[[el]], track.height = 0.08,window.size = 1e7,overlap=FALSE,count_by = "number")
              
      }
            
    }
          
    #ADD BNDs
    
    if (info$size != 0) {

      for (i in 1:nrow(BNDs)) {
        circos.link(BNDs[i,]$chr,as.numeric(BNDs[i,]$pos1),BNDs[i,]$chr2,as.numeric(BNDs[i,]$pos2))
      }

    }
        
    circos.clear()  
  }
} else {

  tool<-args[2]
  outdir<-file.path(args[3])
  pdf(file.path(outdir, paste0("SI00001.", tool, ".svs.pdf")), height=6, width=14)

  for (dir in c(args[1])) {


    bedpe<-file.path(dir, "SI00001.bnd.bedpe") #this can exists or not
    info<-file.info(bedpe)

    if (info$size != 0 && ! is.na(info$size)) { #this is not always true

      BNDs<-fread(bedpe, sep="\t", header=FALSE)
      colnames(BNDs)<-c("chr", "pos1", "end1", "chr2", "pos2", "end2", "type", "num", "s1", "s2")

    }
          
    bed<-file.path(dir, "SI00001.nobnd.bed") #this is always true
    NOBNDs<-fread(bed, sep="\t", header=FALSE)
    colnames(NOBNDs)<-c("chr", "pos", "end", "type")
    NOBNDRanges<-makeGRangesFromDataFrame(NOBNDs, keep.extra.columns = TRUE,start.field = "pos", end.field = "end")
    chrSizes<-seqlengths(Hsapiens)[c(1:24)]
    windows<-tileGenome(chrSizes, tilewidth=1e7, cut.last.tile.in.chrom=T)
    svtypes<-c("INS", "DEL", "DUP", "INV")
    finvalues<-list()
          
    for (l in 1:length(windows)) {
            
      w<-windows[l]
      chr<-as.character(w@seqnames@values)
      start<-as.numeric(w@ranges@start)
      end<-as.numeric(data.frame(w@ranges)$end)
      bps<-as.numeric(w@ranges@width)
      res<-findOverlaps(w, NOBNDRanges)
      idx<-res@to
      subw<-NOBNDs[idx,]
            
      #iterate over types
            
      values<-list()
            
      for (i in 1:length(svtypes)) {
              
        idxtype<-which(subw$type == svtypes[i])
              
        if (length(idxtype) == 0) {
                
          ratio<-0 #no alt
                
        } else {
                
          subsubw<-subw[idxtype,]      
          diff<-subsubw$end-subsubw$pos
          cums<-cumsum(diff)[length(diff)]
          ratio<-cums/bps
                
        }
              
        values[[svtypes[i]]]<-ratio
              
      }
            
      block<-do.call(c, values)
      vector<-data.frame('chr' = chr, 'start' = start, 'end' = end, 'ins' = as.numeric(block[1]), 'del' = as.numeric(block[2]), 'dup' = as.numeric(block[3]), 'inv' = as.numeric(block[4]))
      finvalues[[l]]<-vector
            
    }
          
    heattable<-do.call(rbind,finvalues)
    heatdataframe<-split(data.table(heattable), by="chr")
    reslist<-list()

    for (el in names(heatdataframe)) {
          
      if (min(heatdataframe[[el]]$ins) == max(heatdataframe[[el]]$ins)) {
              
        inscolors<-rep("#FFFFFF",nrow(heatdataframe[[el]]))
              
      } else {
            
        ii<-cut(heatdataframe[[el]]$ins, breaks = seq(min(heatdataframe[[el]]$ins), max(heatdataframe[[el]]$ins), len = 100), include.lowest = TRUE)
        inscolors <- colorRampPalette(c("white", "darkblue"))(99)[ii]
              
      }

      if (min(heatdataframe[[el]]$del) == max(heatdataframe[[el]]$del)) {
              
        delcolors<-rep("#FFFFFF",nrow(heatdataframe[[el]]))
              
      } else {
              
        ii<-cut(heatdataframe[[el]]$del, breaks = seq(min(heatdataframe[[el]]$del), max(heatdataframe[[el]]$del), len = 100), include.lowest = TRUE)
        delcolors <- colorRampPalette(c("white", "mediumorchid4"))(99)[ii]
              
      }
            
      if (min(heatdataframe[[el]]$dup) == max(heatdataframe[[el]]$dup)) {
              
        dupcolors<-rep("#FFFFFF",nrow(heatdataframe[[el]]))
              
      } else {
              
        ii<-cut(heatdataframe[[el]]$dup, breaks = seq(min(heatdataframe[[el]]$dup), max(heatdataframe[[el]]$dup), len = 100), include.lowest = TRUE)
        dupcolors <- colorRampPalette(c("white", "gold1"))(99)[ii]
              
      }
            
      if (min(heatdataframe[[el]]$inv) == max(heatdataframe[[el]]$inv)) {
              
        invcolors<-rep("#FFFFFF",nrow(heatdataframe[[el]]))
              
      } else {
              
        ii<-cut(heatdataframe[[el]]$inv, breaks = seq(min(heatdataframe[[el]]$inv), max(heatdataframe[[el]]$inv), len = 100), include.lowest = TRUE)
        invcolors <- colorRampPalette(c("white", "darkgreen"))(99)[ii]
              
      }
            
      df<-cbind(inscolors,delcolors,dupcolors,invcolors)
      reslist[[el]]<-df
    
    }
         
    mat_col<-matrix(do.call(rbind,reslist),ncol = 4)
          
    circos.par("start.degree" = 90)
    circos.initializeWithIdeogram(species = "hg38")
    circos.genomicHeatmap(heattable, col = mat_col,line_col = as.numeric(factor(heatdataframe[[1]])),connection_height = NULL)
          
    NONBNSs_colors<-list()
    NONBNSs_colors[["INS"]]<-"#39558CFF"
    NONBNSs_colors[["DEL"]]<-"#440154FF"
    NONBNSs_colors[["DUP"]]<-"#FDE725FF"
    NONBNSs_colors[["INV"]]<-"#74D055FF"
          
    #PLOT NON-BNDs first
          
    bed_list<-split(NOBNDs,by="type")
          
    for (el in svtypes) {
            
      if (el %in% names(bed_list)) {
              
        circos.genomicDensity(bed_list[[el]], col = NONBNSs_colors[[el]], track.height = 0.08,window.size = 1e7,overlap=FALSE,count_by = "number")
              
      }
            
    }
          
    #ADD BNDs
    
    if (info$size != 0 && ! is.na(info$size)) {

      for (i in 1:nrow(BNDs)) {
        circos.link(BNDs[i,]$chr,as.numeric(BNDs[i,]$pos1),BNDs[i,]$chr2,as.numeric(BNDs[i,]$pos2))
      }

    }
        
    circos.clear()  
  }  

}
  
#title(main=paste0("SVs - ", toupper(tool)), line = -2, outer = TRUE)
draw(lgd_list_horizontal, y = unit(1, "npc") - circle_size, just = "bottom")

dev.off()


if (file.exists("Rplots.pdf")) {

  file.remove("Rplots.pdf")

}



