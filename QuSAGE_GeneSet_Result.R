## Copyright (C) 2019-2020 Novelbio
## Running Platform: Ubuntu / Novelbrain
## Updator: LI, Yadang          Update time: 2020.01.10
## Author:  LI, Yadang          Create time: 2019.11.27
## Instruction: The function of this script is to create GeneSets results(Data table / Density Curve plot / Confidence Interval plot) based on qs.rds file with QuSAGE algorithm result

## Function list:
##   LoadInputPara
##   get_list_gmt
##   get_cset
##   SetPlotArg
##   plotDC <- plotDensityCurve
##   plotCI <- plotCIs
##   pdfScaleFactor <- pdfScaleFactor
##   GetResultFiles
##     CreateGeneSetInfoTable
##     CreateGeneSetDCplot
##     CreateGeneSetCIplot
##   SetResultFileName
##   ProcessQSresults
##   main

## Global Variable list:
##   crtwd (current work directory)
##   file_gmt (geneset file)
##   curve.num
##   font.size
##   eset
##   ccls
##   cset
##   genesets
##   gset
##   color01
##   color02
##   color10
##   color80
##   label_gmt
##   table_name
##   plotDC_name
##   plotCI_name

library(qusage)
library(getopt)
library(optparse)

## Func: Get input arguments
LoadInputPara <- function () {
  command = matrix(c('outpath'       , 'o', 1, "character",
                     'outinpath'     , 'i', 1, "character",
                     'binpath'       , 'b', 2, "character",
                     'gsetpath'      , 'g', 2, "character",
                     'input_table'   , 't', 2, "character",
                     'input_MSigDB'  , 'm', 2, "character",
                     'input_Novelbio', 'n', 2, "character",
                     'curvenum'      , 'c', 2, "integer",
                     'fontsize'      , 'f', 2, "double"
                    ), byrow = TRUE, ncol = 4)
  args = getopt(command)
  curve.num <<- 50; if (!is.null(args$curvenum)) curve.num <<- args$curvenum
  font.size <<- 3;  if (!is.null(args$fontsize)) font.size <<- args$fontsize
  out.path  <<- args$outpath; print(paste("OutPath:", out.path, sep = " ")); setwd(out.path); crtwd <<- getwd();
  task.path <<- args$outinpath; print(paste("TaskPath:", task.path, sep = " "))
  gset.path <<- args$gsetpath; print(paste("GsetPath:", gset.path, sep = " "))
}

## Function: get genesets list from tmp/gslist.txt
get_list_gmt <- function () {
  file <- paste0(task.path, "tmp/gslist.txt"); #print(file)
  list_gmt <<- as.matrix(read.table(file, he = F)); #print(list_gmt)
}

## Function: get cluster list from tmp/clist.txt
get_cset <- function () {
  file <- paste0(task.path, "tmp/clist.txt"); #print(file)
  cset <<- as.matrix(read.table(file, he = F)); #print(cset)
}

## Function: Create color schemes for plot
SetPlotArg <- function () {
  color01 <<- rainbow(90)
  color02 <<- c("black", rainbow(60))
  color80 <<- c("#00008B", "#FFA500", "#00CD00", "#FF0000", "#9B30FF", "#8B4513",
     "#FFB5C5", "#A9A9A9", "#FFFF00", "#00BFFF", "#006400", "#0000FF", "#8B008B",
     "#8B0000", "#8B8B00", "#76EEC6", "#FF00FF", "#FFFACD", "#363636", "#FF6347",
     "#8B3A62", "#FF1493", "#00688B", "#708090", "#ADFF2F", "#9370DB", "#FF7F00",
     "#CD69C9", "#FFDAB9", "#90EE90", "#6495ED", "#CD3333", "#D2B48C", "#68228B",
     "#FF69B4", "#0000CD", "#458B00", "#FFF68F", "#00FFFF", "#8B4C39", "#1E90FF",
     "#F4A460", "#CAFF70", "#CDCDB4", "#EE6363", "#E066FF", "#EED5D2", "#CD3278",
     "#CDCD00", "#009ACD", "#FF4500", "#228B22", "#BDB76B", "#FFFFE0", "#D2691E",
     "#CDAD00", "#98F5FF", "#00EE76", "#8968CD", "#DCDCDC", "#8B8970", "#FF1493",
     "#0000FF", "#8B7B8B", "#7FFF00", "#D8BFD8", "#EE9A00", "#EE6A50", "#E6E6FA",
     "#008B8B", "#9932CC", "#00FA9A", "#B0C4DE", "#FF6A6A", "#FFD700", "#458B00",
     "#4876FF", "#8A2BE2", "#EEDD82", "#4682B4")
##               "darkblue",        "orange",         "green3",         "red",           "purple",        "saddlebrown",
## "pink1",      "darkgrey",        "yellow",         "deepskyblue",    "darkgreen",     "blue",          "darkmagenta",
## "red4",       "yellow4",         "aquamarine2",    "magenta1",       "lemonchiffon1", "grey21",        "tomato1",
## "hotpink4",   "deeppink1",       "deepskyblue4",   "slategrey",      "greenyellow",   "mediumpurple",  "darkorange1",
## "orchid3",    "peachpuff1",      "lightgreen",     "cornflowerblue", "brown3",        "tan",           "darkorchid4",
## "hotpink",    "blue3",           "chartreuse4",    "khaki1",         "cyan",          "salmon4",       "dodgerblue1",
## "sandybrown", "darkolivegreen1", "lightyellow3",   "indianred2",     "mediumorchid1", "mistyrose2",    "violetred3",
## "yellow3",    "deepskyblue3",    "orangered1",     "forestgreen",    "darkkhaki",     "lightyellow",   "chocolate",
## "gold3",      "cadetblue1",      "springgreen2",   "mediumpurple3",  "gainsboro",     "lemonchiffon4", "deeppink",
## "blue1",      "thistle4",        "chartreuse",     "thistle",        "orange2",       "coral2",        "lavender",
## "cyan4",      "darkorchid",      "medspringgreen", "lightsteelblue", "indianred1",    "gold",          "chartreuse4",
## "royalblue1", "blueviolet",      "lightgoldenrod", "steelblue"
  color10 <<- c(color80[1:10], rep("grey", 70))
}

## Instruction: This function is modified code based on plotDensityCurve from qusage library
plotDC <- function (QSarray, path.index = 1:numPathways(QSarray), zeroLine = TRUE,
		    addVIF = !is.null(QSarray$vif), col = NULL, plot = TRUE,
		    add = FALSE, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL,
		    type = "l", ...) {
    if (is.character(path.index)) {
        path.index = match(path.index, names(QSarray$pathways))
    }
    if (all(is.na(QSarray$path.mean[path.index]))) {
        stop("no non-missing pathways in path.index")
    }
    scaleFactor = pdfScaleFactor(QSarray, addVIF = addVIF)
    if (is.null(xlim)) {
        xlim = range(calcBayesCI(QSarray, addVIF = addVIF)[, path.index], na.rm = T)
        xlim = xlim + (xlim - mean(xlim)) * 2
    }
    if (is.null(ylim)) {
        ylim = range(t(QSarray$path.PDF[, path.index])/scaleFactor[path.index], na.rm = T)
    }
    if (is.null(xlab))  xlab = "Gene Set Activity"
    if (is.null(ylab))  ylab = "Density"
    if (is.null(col))   col = par("col")
    if (!add & plot) {
        plot(0, type = "n", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, ...)
    }
    if (length(col) != length(path.index)) {
        col = rep(col, length.out = length(path.index))
    }
    if (zeroLine & plot) abline(v = 0, lty = 2)
    retVal = list()
    for (i in 1:length(path.index)) {
        path = path.index[i]
        x = getXcoords(QSarray, path, addVIF = addVIF)
        y = QSarray$path.PDF[, path]/scaleFactor[path]
        if (plot) {
            lines(x, y, col = col[i], type = type, ...)
        }
        retVal[[i]] = data.frame(x, y)
    }
    names(retVal) = colnames(QSarray$path.PDF)[path.index]
    invisible(retVal)
}

## Instruction: This function is modified code based on plotCIs from qusage library
plotCI <- function (QSarray, path.index = 1:numPathways(QSarray), sort.by = c("mean", "p", "none"),
    lowerBound = 0.025, upperBound = 1 - lowerBound,
    col = NULL, use.p.colors = TRUE, p.breaks = NULL, p.adjust.method = "fdr",
    addLegend = use.p.colors, lowerColorBar = "none", lowerColorBar.cols = NULL,
    addGrid = TRUE, x.labels = NULL, cex.xaxis = 1, shift = 0,
    add = FALSE, ylim = NULL, xlim = NULL, ylab = NULL, xlab = NULL,
    main = NULL, sub = NULL, type = "p", ...)
{
    if (is.character(path.index)) path.index = match(path.index, names(QSarray$pathways))
    if (all(is.na(QSarray$path.mean[path.index]))) stop("no non-missing pathways in path.index")
    means = QSarray$path.mean[path.index]; #print(means)
    CIs = calcBayesCI(QSarray, low = lowerBound, up = upperBound)[, path.index, drop = FALSE]
    p.vals = pdf.pVal(QSarray, direction = F)[path.index]
    p.vals = p.adjust(p.vals, p.adjust.method)
    p.vals = p.vals * sign(means)
    if (is.null(x.labels)) x.labels = names(means); #print(x.labels)
    sort.by = match.arg(sort.by)
    if (sort.by == "mean") { #print(1)
        ord = order(means, decreasing = T, na.last = NA)
    }
    else if (sort.by == "p") { #print(2)
        ord = order(-log10(abs(p.vals)) * sign(p.vals), decreasing = T,
            na.last = NA)
    }
    else { #print(3)
        ord = 1:length(means)
    }; #print(x.labels[ord]); print(ord)
    originalPar = par(no.readonly = T)
    additional_args <- list(...)
    if (!"srt" %in% names(additional_args))
        par(srt = 60)
    par(...)
    if (!is.null(col)) {
        if (length(means)%%length(col) != 0) {
            warning("length of col is not a multiple of input length")
        }
        col = rep(col, length.out = length(means))
    }
    else {
        col = rep(par("col"), length.out = length(means))
    }
    if (is.null(p.breaks)) {
        p.breaks = c(0.001, 0.005, 0.01, 0.05, 0.1)
    }
    else {
        p.breaks = abs(p.breaks)
        p.breaks = as.numeric(names(table(p.breaks)))
        p.breaks = p.breaks[p.breaks < 1 & p.breaks > 0]
        p.breaks = p.breaks[order(p.breaks)]
    }
    br.ln = length(p.breaks)
    p.breaks.twoSided <- c(-1, -rev(p.breaks), 0, p.breaks, 1)
    if (use.p.colors) {
        p.colorScheme <- c(rgb(0, seq(0, 1, length.out = br.ln + 1), 0), rgb(seq(1, 0, length.out = br.ln + 1), 0, 0))
        bar.col = p.colorScheme[findInterval(p.vals, p.breaks.twoSided, rightmost.closed = T)]
        bar.col[which(p.vals == 0)] = c("#00FF00", "#FF0000")[(means > 0) + 1][which(p.vals == 0)]
    }
    else {
        bar.col = col
    }
    if (is.null(xlab)) xlab = ""
    if (is.null(ylab)) ylab = "Gene Set Activity"
    if (!add) {
        if (is.null(ylim)) {
            ylim = range(CIs, na.rm = T); ylim[2] = ylim[2] + (ylim[2] - ylim[1]) * 0.125 ## set ylim to change the position of colour bar
        }
        plot(means[ord], type = "n", las = 1, ylim = ylim, axes = FALSE, xlab = xlab, ylab = ylab, main = main, ...)
        if (addGrid) {
            abline(v = 1:length(ord), col = gray(seq(0.5, 1, length.out = 10)), lty = 2)
        }
        axis(2, las = 1, ...)
        if (is.null(list(...)$xaxt) || list(...)$xaxt != "n") {
            Oldcex = par("cex")
            par(cex = 0.5 * Oldcex * cex.xaxis)
            axis(1, at = ord, las = 2, labels = rep("", length(ord)),  ...)
            Ys <- par("usr")[3] - par()$cxy[2] * par()$mgp[2]
            text(1:length(ord), Ys, adj = 1, labels = x.labels[ord], xpd = TRUE, ...); #print(x.labels[ord]); print(ord)
            par(cex = Oldcex)
        }
        abline(h = 0, lty = 2)
    }
    arrows(1:length(ord) + shift, CIs[2, ord], 1:length(ord) + shift, CIs[1, ord], code = 3, length = 0.1, angle = 90, col = bar.col[ord])
    points(1:length(ord) + shift, means[ord], type = type, col = col[ord],  ...)
    if (is.null(x.labels))
        x.labels = names(QSarray$path.mean)[path.index]
    if (use.p.colors && addLegend && !add) {
        p.colorScheme <- c(rgb(0, seq(1, 0, length.out = br.ln + 1), 0), rgb(seq(0, 1, length.out = br.ln + 1), 0, 0))
        p.labels <- round(c(-p.breaks, 1, rev(p.breaks)), 3)
        n = length(p.labels)
        usr = par("usr"); #print(usr)
        strsize = strwidth("W")
        yvalues <- usr[4] - c(0, strheight("W") * 1.5); #print(yvalues)
        xvalues <- usr[2] - (seq(strsize * (n + 2), 0, length.out = n + 2) * 1.2); #print(xvalues)
        text(mean(xvalues), yvalues[0] + strheight("X"), labels = "P-values", pos = 3, xpd = T)
        for (j in 1:(n + 1)) {
            polygon(c(xvalues[j], xvalues[j + 1], xvalues[j + 1], xvalues[j]), c(yvalues[1], yvalues[1], yvalues[2], yvalues[2]), col = p.colorScheme[j], border = p.colorScheme[j])
        }
        for (j in 1:n) {
            text(xvalues[j + 1], yvalues[2] - strheight("W"), p.labels[j], adj = 1, srt = 90)
        }
    }
    lowBar.vals = NULL
    if (!is.character(lowerColorBar)) {
        if (length(lowerColorBar) != length(path.index)) {
            warning("lowerColorBar not the same length as path.index. Color bar omitted")
        }
        else {
            lowBar.vals = lowerColorBar
            if (is.null(lowerColorBar.cols)) {
                lowerColorBar.cols = rgb(seq(1, 0, length.out = 6),
                  seq(1, 0, length.out = 6), 1)
            }
            lowBar.breaks = seq(range(lowBar.vals)[1], range(lowBar.vals)[2],
                length.out = length(lowerColorBar.cols) + 1)
        }
    }
    else if (lowerColorBar == "absolute") {
        if (is.null(QSarray$absolute.p)) {
            warning("Absolute p-values not found. Color bar omitted.")
        }
        else {
            lowBar.vals = QSarray$absolute.p[path.index]
            lowBar.vals <- p.adjust(abs(lowBar.vals), method = p.adjust.method)
            lowBar.breaks = p.breaks.twoSided
            if (is.null(lowerColorBar.cols)) {
                lowerColorBar.cols = c(rgb(0, seq(0, 1, length.out = br.ln + 1), 0, 0), rgb(seq(1, 0, length.out = br.ln + 1), 0, 0))
            }
        }
    }
    else if (lowerColorBar == "homogeneity") {
        if (is.null(QSarray$homogeneity)) {
            warning("Homogeneity score not found. Color bar omitted.")
        }
        else {
            lowBar.vals = QSarray$homogeneity[path.index]
            if (is.null(lowerColorBar.cols)) {
                lowerColorBar.cols = rgb(seq(1, 0, length.out = 6), seq(1, 0, length.out = 6), 1)
            }
            lowBar.breaks = seq(range(lowBar.vals)[1], range(lowBar.vals)[2],
                length.out = length(lowerColorBar.cols) + 1)
        }
    }
    if (!is.null(lowBar.vals)) {
        col = lowerColorBar.cols[findInterval(lowBar.vals, lowBar.breaks, rightmost.closed = T)]
        for (i in 1:length(ord)) {
            usr = par()$usr
            polygon(i + c(-0.5, -0.5, 0.5, 0.5), usr[3] + c(0.01,
                0.03, 0.03, 0.01) * (usr[4] - usr[3]), col = col[ord][i],
                border = col[ord][i])
        }
    }
    box(...)
    par(originalPar[c("srt", "pch", names(list(...)))])
}

## Instruction: This function is the original code from qusage library
pdfScaleFactor <- function(QSarray, addVIF=!is.null(QSarray$vif)){
  if(is.null(QSarray$vif) && addVIF){stop("vif is undefined for QSarray object. addVIF can not be set to true.")}
  sif = sapply(1:numPathways(QSarray),function(i){ifelse(addVIF,sqrt(QSarray$vif[i]),1)})
  sif[is.na(sif)] = 1
  pdfSum = colSums(QSarray$path.PDF)

  ##get the pdf range
  ranges = QSarray$ranges * 2

  ##the scale factor is essentially the distance between points in the x coordinates times the (current) sum of the pdf
  scaleFactor = (ranges*sif) / (QSarray$n.points-1) * pdfSum
  scaleFactor
}

## Function: Set output file names inclue DensityCurvesPlot, ConfidenceIntervalPlot and DataTable
## Input:    cluster (cluster name)
SetResultFileName <- function (cluster) {
  label <- paste("Cluster", cluster, label_gmt, sep = "_")
  outwd <- paste0(gset.path, "/", label_gmt)
  if (!file.exists(outwd)) dir.create(outwd, recursive = T)
  table_name  <<- paste0(outwd, "/", "GeneSetInfo_", label, ".txt")
  plotDC_name <<- paste0(outwd, "/", "GeneSetDCplot_", label, ".png")
  plotCI_name <<- paste0(outwd, "/", "GeneSetCIplot_", label, ".png")
}

## Function: Create qs.result and put it into a RDS file
## Input:    cluster (cluster name)
ProcessQSresults <- function (cluster) { #print(paste("cluster:",cluster))
  file_qs_result <- paste0(task.path, "tmp/", label_gmt, "/Cluster_", cluster, ".qs.rds");
  if (!file.exists(file_qs_result)) print(paste("Cannot find file:", file_qs_result))
  qs.result <<- readRDS(file_qs_result)
}

## Function: Create DataTable based on qs.result
CreateGeneSetInfoTable <- function (
    table_name,  ## file name of DataTable
    table_info   ## info of DataTable
  ) {
  colnames(table_info)[1] <- "GeneSet.name"
  write.table(table_info, file = table_name, sep = "\t", row.names = F)
  print(paste("Output:", table_name))
}

## Function: Create DensityCurvesPlot based on qs.result
CreateGeneSetDCplot <- function (
    plotDC_name, ## file name of DensityCurvesPlot
    qs.result,   ## output of qusage
    plot.order,  ## order of items will be plot
    top.num,     ## nummer of curves will be plot
    clusID,      ## cluster name
    genesets     ## names of genesets
  ) {
## plotDC legend order setting
  ordered.legend = character()
  for (i in 1:length(plot.order)) {
    #ordered.legend[i] <- names(genesets)[(plot.order[i])]
    ordered.legend[i] <- genesets[(plot.order[i])]
    if (i > 0) next
    #print(paste(plot.order[i],names(genesets)[plot.order[i]],ordered.legend[i],sep=" "))
    print(paste(plot.order[i],genesets[plot.order[i]],ordered.legend[i],sep=" "))
  }
  #print(names(genesets))
  #print(ordered.legend)

## PlotDC setting
  png(file=plotDC_name, width = 1800, height = 1000)
    par(mgp = c(3,1,0), mar=c(6,6,3,2))
    plotDC(qs.result, path.index=plot.order[top.num:1], lwd=3, col=color80[top.num:1],
           cex.axis=2, cex.lab=2, main=clusID, cex.main=2)
    #plotDC(qs.result, path.index=plot.order[1], lwd=9, col=color80[1],
    #       cex.axis=2, cex.lab=2,
    #       add=T)
    legend("topleft", legend = ordered.legend[1:top.num], text.col = color80, bty="n")#xpd=T)
  dev.off()

#  pdf(file=sub("png", "pdf", plotDC_name), width = 27, height = 15)
#  par(mgp = c(3,1,0), mar=c(6,6,3,3))
#  plotDensityCurves(qs.result, path.index=plot.order[top.num:1], lwd=5, col=colors80[top.num:1],
#		    cex.axis=2, cex.xaxis=2, cex.lab=2,
#		    xlab="Gene Set Activity", ylab="Density", main=clusID, cex.main=2)
#  plotDensityCurves(qs.result, path.index=plot.order[1],
#		    lwd=10, col=colors80[1], cex.axis=2, cex.xaxis=2, cex.lab=2,
#		    add=T)
#  legend("topleft", legend = ordered.legend[1:top.num], text.col = colors80,
#	 bty="o", bg="transparent") #xpd = T)
#  dev.off()

}

## Function: Create ConfidenceIntervalPlot based on qs.result
CreateGeneSetCIplot <- function (
    plotCI_name, ## file name of ConfidenceIntervaPlot
    qs.result,   ## output of qusage
    plot.order,  ## order of items will be plot
    top.num,     ## nummer of curves will be plot
    clusID,      ## cluster name
    genesets     ## names of genesets
  ) {
  max.char <- 30;

## PlotCI arguments setting
  p.vals <- pdf.pVal(qs.result); class(p.vals); #print(p.vals)
  topSets <- order(p.vals); #print(topSets)#print(topSets[1:25])
  path.index = topSets[1:top.num]
  means = qs.result$path.mean[path.index]
  ord = order(means, decreasing = T, na.last = NA)
  x.labels = names(means); #print(x.labels); #for (i in 1:length(ord)) print(c(ord[i], x.labels[ord[i]]))
  labelength = nchar(x.labels); #print(labelength); print(max(labelength))
  max.length = max(labelength)
  if (max.length > max.char) max.length <- max.char; #print(max.length)
  qs.CIplot <- qs.result
  for (n in 1:length(names(qs.CIplot$path.mean))){
    if (nchar(names(qs.CIplot$path.mean)[n]) > max.char) {
      lab27 <- substr(names(qs.CIplot$path.mean)[n], 1, 27)
      lab30 <- paste0(lab27, "..."); #print(lab30)
      names(qs.CIplot$path.mean)[n] <- lab30; #print(names(qs.CIplot$path.mean)[n])
    }
  };  #print(names(qs.CIplot$path.mean))
  #print(path.index)
  #print(plot.order[1:top.num])
  mar.d = 1 + max.length * font.size * 0.24;         #print(paste('mar.d:', mar.d))

  #mar.l = 5 + labelength[ord[1]] * font.size * 0.15;
  mar.l = 5 + max.length * font.size * 0.15;
  if (mar.l < font.size * 3) mar.l = font.size * 3;       #print(paste('mar.l:', mar.l))

  mar.u = font.size + 1;                                  #print(paste('mar.u:', mar.u))

  plotCI.width = 1800;
  plotCI.width = 300 + 30 * top.num + 10 * mar.l
  if (font.size > 4) plotCI.width = 300 + 6.2 * top.num * font.size + 10 * mar.l; #print(paste('width:', plotCI.width))

  plotCI.height = 1000;
  plotCI.height = 1000 + mar.l * 10;                      #print(paste('height:', plotCI.height))

  mgp1 = font.size + 3;                                   #print(paste('mgp1:', mgp1))

## PlotCI setting
  png(file=plotCI_name, width = plotCI.width, height = plotCI.height)
    par(mgp = c(mgp1,1,0))
    plotCI(qs.CIplot, path.index=plot.order[1:top.num], mar=c(mar.d, mar.l, mar.u, 1), lwd=3, #col=colors2[1:top.num],
           cex.axis=0.5 * font.size, cex.xaxis=font.size, cex.lab= 0.75 * font.size,
           main=clusID, cex.main= 0.75 * font.size)
  dev.off()

}

## Function: Create QuSAGE results include DensityCurvesPlots, ConfidenceIntervalPlots and DataTable based on qs.result
GetResultFiles <- function (
    qs.result,   ## output of qusage
    cluster,     ## cluster ID
    table_name,  ## file name of DataTable
    plotDC_name, ## file name of DensityCurvesPlot
    plotCI_name  ## file name of ConfidenceIntervaPlot
  ) {
  if (!is.numeric(cluster)) clusID <- cluster
  if ( is.numeric(cluster)) clusID <- paste("Cluster", cluster)

  genesets <- names(qs.result$pathway)

## plotDC and plotCI curve num setting
  top.num = length(genesets);
  if ((curve.num != 0) & (top.num > curve.num)) top.num <- curve.num

## output genesets info table
## Change the order of gene sets in QuSAGE result table from large logfoldchange to small in the condition of P-value equals to zero.
  table <- qsTable(qs.result, number = (length(genesets))); #print(table) #print(dim(table)) ## class(table): data.frame ## original table
  tableA <- table[!is.na(table$log.fold.change),]; #print(tableA); #print(dim(tableA))
  #tableB <- order(tableA$p.Value, -abs(tableA$log.fold.change))
  tableU <- tableA[tableA$p.Value == 0,]; #print(dim(tableU))
  tableD <- tableA[tableA$p.Value != 0,]; #print(dim(tableD))
  tableO <- tableU[order(abs(tableU$log.fold.change), decreasing = T), ]
  tableB <- rbind(tableO, tableD); #print(tableB); #print(dim(tableB))
  if (is.na(tableO$p.Value[1])) tableB <- tableD;
  tableT <- tableB[1:top.num,]; #print(tableT); #print(dim(tableT))
  tableC <- tableT[order(tableT$log.fold.change, decreasing = T), ]; #print(tableC) #print(dim(tableC));
  plot.order <- as.numeric(row.names(tableC)); #print("plot.order:"); print(plot.order)

  table0 <- tableU[order(tableU$log.fold.change, decreasing = T), ]
  tableZ <- rbind(table0, tableD)
  if (is.na(table0$p.Value[1])) tableZ <- tableD;
  CreateGeneSetInfoTable(table_name, tableZ)
  CreateGeneSetDCplot(plotDC_name, qs.result, plot.order, top.num, clusID, genesets)
  CreateGeneSetCIplot(plotCI_name, qs.result, plot.order, top.num, clusID, genesets)
}

main <- function () {
  print(date())
  print("Running QuSAGE_GeneSet_Result_20200110.R")
  originalPar = par(no.readonly = T) ## store the original full list of parameters
  task.color <<- 0
  LoadInputPara()
  get_list_gmt()
  get_cset()
  SetPlotArg()
  for (f in list_gmt) {
    label_gmt <<- f; print(label_gmt)
    for (i in unique(cset)) {
      print(paste("Cluster", i, date()))
      SetResultFileName(i)
      ProcessQSresults(i)
      GetResultFiles(qs.result, i, table_name, plotDC_name, plotCI_name)
    }
    date()
  }
}

warnings()
main()
#task.color <- 1 ## code to stop the task for debugging
quit(status = task.color)
