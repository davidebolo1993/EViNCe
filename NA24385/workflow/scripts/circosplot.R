#!/usr/bin/env Rscript

library(circlize)
library(data.table)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 3) { #then, we must plot as a single plot
  
  tool<-args[2]
  outdir<-file.path(args[3])
  pdf(file.path(outdir, paste0("GM24385.", tool, ".svs.pdf")), height=10, width=13)

  if (tool != "truth") {
  
    bedpe<-file.path(args[1], "GM24385.bnd.bedpe")
    BNDs<-fread(bedpe, sep="\t", header=FALSE)
    BNDs$V1<-paste0("chr", BNDs$V1)
    BNDs$V4<-paste0("chr", BNDs$V4)
    colnames(BNDs)<-c("chr", "pos1", "end1", "chr2", "pos2", "end2", "type", "num", "s1", "s2")
    
    bed<-file.path(args[1], "GM24385.nobnd.bed")
    NOBNDs<-fread(bed, sep="\t", header=FALSE)
    NOBNDs$V1<-paste0("chr", NOBNDs$V1)
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
        subsubw<-subw[idxtype,]
        
        if (nrow(subsubw) == 0) {
          
          ratio<-0
          
        } else {
          
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
    circos.initializeWithIdeogram(species = "hg19")
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
    
    for (i in 1:nrow(BNDs)) {
      circos.link(BNDs[i,]$chr,as.numeric(BNDs[i,]$pos1),BNDs[i,]$chr2,as.numeric(BNDs[i,]$pos2))
    }
    
    circos.clear()
    
    title(main=paste0("SVs - ", tool))
    legend(1.2,0.5, names(NONBNSs_colors),col=as.character(sapply(NONBNSs_colors, function(x) x[1])),pch=c(16,16,16,16),cex=0.75,title="SVTYPE",bty='n')
    legend(1.2,0.3,legend="TRA",col="black",lty=1,cex=0.75,lwd=1.2,bty='n')

  } else {

    bed<-file.path(args[1], "GM24385.nobnd.bed")
    NOBNDs<-fread(bed, sep="\t", header=FALSE)
    NOBNDs$V1<-paste0("chr", NOBNDs$V1)
    colnames(NOBNDs)<-c("chr", "pos", "end", "type")
    NOBNDRanges<-makeGRangesFromDataFrame(NOBNDs, keep.extra.columns = TRUE,start.field = "pos", end.field = "end")
    chrSizes<-seqlengths(Hsapiens)[c(1:24)]
    windows<-tileGenome(chrSizes, tilewidth=1e7, cut.last.tile.in.chrom=T)
    svtypes<-c("INS", "DEL")
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
        subsubw<-subw[idxtype,]
        
        if (nrow(subsubw) == 0) {
          
          ratio<-0
          
        } else {
          
          diff<-subsubw$end-subsubw$pos
          cums<-cumsum(diff)[length(diff)]
          ratio<-cums/bps
          
        }
        
        values[[svtypes[i]]]<-ratio
        
      }
      
      block<-do.call(c, values)
      vector<-data.frame('chr' = chr, 'start' = start, 'end' = end, 'ins' = as.numeric(block[1]), 'del' = as.numeric(block[2]))
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
      
      df<-cbind(inscolors,delcolors)
      reslist[[el]]<-df
    }
    
    mat_col<-matrix(do.call(rbind,reslist),ncol = 2)
    
    circos.par("start.degree" = 90)
    circos.initializeWithIdeogram(species = "hg19")
    circos.genomicHeatmap(heattable, col = mat_col,line_col = as.numeric(factor(heatdataframe[[1]])),connection_height = NULL)
    
    NONBNSs_colors<-list()
    NONBNSs_colors[["INS"]]<-"#39558CFF"
    NONBNSs_colors[["DEL"]]<-"#440154FF"
    
    #PLOT NON-BNDs first
    
    bed_list<-split(NOBNDs,by="type")
    
    for (el in svtypes) {
      
      if (el %in% names(bed_list)) {
        
        circos.genomicDensity(bed_list[[el]], col = NONBNSs_colors[[el]], track.height = 0.08,window.size = 1e7,overlap=FALSE,count_by = "number")
        
      }
      
    }
        
    circos.clear()
    
    title(main=paste0("SVs - ", tool))
    legend(1.2,0.5, names(NONBNSs_colors),col=as.character(sapply(NONBNSs_colors, function(x) x[1])),pch=c(16,16,16,16),cex=0.75,title="SVTYPE",bty='n')

  }

  dev.off()

} else {
  
  tool<-args[3]
  outdir<-file.path(args[4])
  pdf(file.path(outdir, paste0("GM24385.", tool, ".svs.pdf")), height=10, width=15)
  layout(matrix(1:3, 1,3), widths=c(3,1,3))

  for (dir in c(args[1], "empty", args[2])) {

    if (dir == "empty") {

      plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
      legend(0,0.6,names(NONBNSs_colors),col=as.character(sapply(NONBNSs_colors, function(x) x[1])),pch=c(16,16,16,16),cex=0.75,title="SVTYPE",bty='n')
      legend(0,0.54,legend="TRA",col="black",lty=1,cex=0.75,lwd=1.2,bty='n', seg.len=0.2)
      next

    }
    
    bedpe<-file.path(dir, "GM24385.bnd.bedpe")
    BNDs<-fread(bedpe, sep="\t", header=FALSE)
    BNDs$V1<-paste0("chr", BNDs$V1)
    BNDs$V4<-paste0("chr", BNDs$V4)
    colnames(BNDs)<-c("chr", "pos1", "end1", "chr2", "pos2", "end2", "type", "num", "s1", "s2")
      
    bed<-file.path(dir, "GM24385.nobnd.bed")
    NOBNDs<-fread(bed, sep="\t", header=FALSE)
    NOBNDs$V1<-paste0("chr", NOBNDs$V1)
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
        subsubw<-subw[idxtype,]
          
        if (nrow(subsubw) == 0) {
            
          ratio<-0
            
        } else {
            
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
    circos.initializeWithIdeogram(species = "hg19")
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
      
    for (i in 1:nrow(BNDs)) {
      circos.link(BNDs[i,]$chr,as.numeric(BNDs[i,]$pos1),BNDs[i,]$chr2,as.numeric(BNDs[i,]$pos2))
    }
    
    circos.clear() 

  }
  
  title(main=paste0("SVs - ", tool), line = -5, outer = TRUE) 
  dev.off()

}


if (file.exists("Rplots.pdf")) {

  file.remove("Rplots.pdf")

}

