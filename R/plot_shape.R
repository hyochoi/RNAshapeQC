#' Plot gene body coverage
#'
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param geneNames gene names per file. If NULL, Gene i with the same length of pileupPath be set. Default is NULL.
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @param scale TRUE/FALSE returns the scaled/unscaled normalized transcript coverage. Default is TRUE.
#' @param stat 1 and 2 return median and mean normalized coverage curves per sample, respectively. Default is 1.
#' @param plot TRUE/FALSE turns on/off the normalized transcript coverage plot. Default is TRUE.
#' @param sampleInfo a sample information table including sample id. The number of rows is equal to the number of samples.
#' @return a matrix and a plot, or a matrix for the gene body coverage where plot is TRUE or FALSE, respectively.
#' @importFrom dplyr mutate select inner_join
#' @importFrom grDevices colorRampPalette
#' @importFrom ggplot2 ggplot aes geom_line labs
#' @importFrom ggplot2 scale_colour_gradientn theme element_rect element_text
#' @export

plot_GBC = function(pileupPath, geneNames, rnum=100, method=1, scale=TRUE, stat=2, plot=TRUE, sampleInfo) {

  scale.log.normlist = scale_pileup.list(pileupPath, geneNames, rnum=rnum, method=method, scale=scale)
  scale.arr <- simplify2array(scale.log.normlist)

  if (stat==1) {
    # Stat 1: Median curve per sample
    scale.geom <- matrix(matrixStats::rowMedians(matrix(scale.arr, nrow=prod(dim(scale.arr)[1:2]), ncol=dim(scale.arr)[3])), nrow=dim(scale.arr)[1])

  } else if (stat==2) {
    # Stat 2: Mean curve per sample
    scale.geom <- apply(scale.arr, 1, base::mean)

  } else {
    stop(stat," is not an option for stat.")
  }

  rownames(scale.geom) <- c(seq_len(rnum))
  colnames(scale.geom) <- colnames(scale.log.normlist[[1]])

  lgd <- sampleInfo %>%
    mutate(RatioIntron=INTRONIC_BASES/CODING_BASES) %>%
    select(SampleID, RatioIntron)

  GBP <- convert_Chr2Numcol(convert_pivot.longer(scale.geom, c("region", "sample", "scale.geom")) %>%
                              mutate(sample=gsub("\\.", "-", sample)) %>%
                              select(region, sample, scale.geom) %>%
                              inner_join(lgd, by=c("sample"="SampleID")), 1)

  p <- NULL
  if (plot) {
    myPalette <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))

    p <- ggplot2::ggplot(GBP, aes(x=region, colour=RatioIntron, group=sample)) +
      geom_line(aes(y=scale.geom), alpha=1, show.legend=TRUE) +
      labs(title="", x="Gene body percentile (5'\u21923')", y=paste0(c("Median", "Mean")[stat]," scaled normalized coverage")) +
      scale_colour_gradientn(colours=myPalette(100), limits=c(min(GBP$RatioIntron), max(GBP$RatioIntron)), name="Ratio intron") +
      theme(legend.position="bottom",
            panel.background=element_rect(fill="gray97"),
            strip.text.x=element_text(size=12, color="black"),
            strip.text.y=element_text(size=12, color="black"),
            strip.background=element_rect(color="NA", fill="white", linewidth=1, linetype="solid"),
            plot.title=element_text(hjust=0.5, face="bold"))
    # print(p)
  }

  return(list(GBP=GBP, plot=p))
}


#' Plot gene body coverage with optimal samples
#'
#' @param stat 1 and 2 return median and mean normalized coverage curves per sample, respectively. Default is 1.
#' @param plot TRUE/FALSE turns on/off the normalized transcript coverage plot. Default is TRUE.
#' @param sampleInfo a sample information table including sample id. The number of rows is equal to the number of samples.
#' @param GBCresult results of the gene body coverage with all samples
#' @param auc.vec a vector with SOI per sample
#' @return a matrix and a plot, or a matrix for the gene body coverage where plot is TRUE or FALSE, respectively.
#' @importFrom dplyr filter select inner_join
#' @export

plot_GBCos = function(stat=2, plot=TRUE, sampleInfo, GBCresult, auc.vec) {

  # Update with optimal samples and PD
  GBP <- GBCresult$GBP %>%
    filter(sample %in% as.vector(auc.vec[auc.vec$SOI=="Optimal", c("Sample")])$Sample) %>%
    inner_join(auc.vec %>% select(Sample, PD), by=c("sample"="Sample"))

  pRI <- NULL
  pPD <- NULL
  if (plot) {
    myPalette <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))

    pRI <- ggplot2::ggplot(GBP, aes(x=region, colour=RatioIntron, group=sample)) +
      geom_line(aes(y=scale.geom), alpha=1, show.legend=TRUE) +
      labs(title="", x="Gene body percentile (5'\u21923')", y=paste0(c("Median", "Mean")[stat]," scaled normalized coverage")) +
      scale_colour_gradientn(colours=myPalette(100), limits=c(min(GBCresult$GBP$RatioIntron), max(GBCresult$GBP$RatioIntron)), name="Ratio intron") + # the original range of legend
      theme(legend.position="bottom",
            panel.background=element_rect(fill="gray97"),
            strip.text.x=element_text(size=12, color="black"),
            strip.text.y=element_text(size=12, color="black"),
            strip.background=element_rect(color="NA", fill="white", linewidth=1, linetype="solid"),
            plot.title=element_text(hjust=0.5, face="bold"))

    pPD <- ggplot2::ggplot(GBP, aes(x=region, colour=PD, group=sample)) +
      geom_line(aes(y=scale.geom), alpha=1, show.legend=TRUE) +
      labs(title="", x="Gene body percentile (5'\u21923')", y=paste0(c("Median", "Mean")[stat]," scaled normalized coverage")) +
      scale_colour_gradientn(colours=myPalette(100), limits=c(min(GBP$PD), max(GBP$PD)), name="PD") +
      theme(legend.position="bottom",
            panel.background=element_rect(fill="gray97"),
            strip.text.x=element_text(size=12, color="black"),
            strip.text.y=element_text(size=12, color="black"),
            strip.background=element_rect(color="NA", fill="white", linewidth=1, linetype="solid"),
            plot.title=element_text(hjust=0.5, face="bold"))

    # print(p)
  }

  return(list(GBP=GBP, plotRI=pRI, plotPD=pPD))
}


#' Plot degraded/intact index outputs
#'
#' @param DR a the number of genes x the number of samples matrix of decay rates
#' @param DIIresult outputs from \code{get_DII} function
#' @param cutoff numeric threshold on projection depth used to classify samples.
#' @param outFile a directory with a file name to save outputs. Default is NULL.
#' @return figures for the distribution of DII by PD; and the heatmap of DR.
#' @references https://jtr13.github.io/cc21fall2/raincloud-plot-101-density-plot-or-boxplotwhy-not-do-both.html
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange mutate select case_when
#' @importFrom ggplot2 ggplot aes position_nudge position_jitter
#' @importFrom ggplot2 geom_point geom_boxplot scale_fill_manual
#' @importFrom ggplot2 scale_color_identity labs guides guide_legend
#' @importFrom ggplot2 theme element_rect element_text coord_flip
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation rowAnnotation
#' @importFrom ComplexHeatmap anno_points Legend draw packLegend
#' @importFrom circlize colorRamp2
#' @importFrom colortools setColors
#' @importFrom ggpubr as_ggplot
#' @importFrom grid grid.grabExpr pushViewport viewport unit gpar
#' @importFrom grDevices png dev.off
#' @importFrom ggplot2 margin
#' @export

plot_DIIwt = function(DR, DIIresult, cutoff=3, outFile=NULL) {

  DR.mat <- DR
  ds.vec <- DIIresult$ds.vec
  gene.df <- DIIresult$gene.df

  ds.vec2 <- ds.vec[, c("DS", "PD")]
  df <- data.frame(var=rep(colnames(ds.vec2), each=nrow(ds.vec2)),
                   value=as.vector(as.matrix(ds.vec2)))
  df$var <- factor(df$var, levels=rev(colnames(ds.vec2)), ordered=TRUE)

  # Distribution of PD by DII
  make_vioplot <- function(df_sub, cutoff) {
    ggplot2::ggplot(df_sub, aes(x=var, y=value, fill=var)) +
      PupillometryR::geom_flat_violin(position=position_nudge(x=0.2), alpha=0.4) +
      geom_point(aes(color=case_when(
        var=="DS" ~ "grey",
        var=="PD" & value>cutoff ~ "#7C3F11",
        var=="PD" & value<=cutoff ~ "mediumaquamarine")),
        position=position_jitter(width=0.15, seed=12345), size=1.5, alpha=0.5, show.legend=F) +
      geom_boxplot(width=0.25, outlier.shape=NA, alpha=0.5) +
      scale_fill_manual(values=c(rep("grey", 2))) +
      scale_color_identity() +
      labs(x="", y="", fill="", title="") +
      guides(fill=guide_legend(nrow=1, byrow=TRUE))+
      theme(legend.position="none",
            legend.margin=margin(-5, 0, 0, 0),
            panel.background=element_rect(fill="gray97"),
            axis.text=element_text(size=25),
            axis.title=ggplot2::element_blank()) +
      coord_flip()
  }

  p1 <- make_vioplot(filter(df, var=="DS"), cutoff)
  p2 <- make_vioplot(filter(df, var=="PD"), cutoff)

  d <- patchwork::wrap_plots(p1, p2, ncol=1)

  # Heatmap of DR
  ds.vec3 <- ds.vec %>%
    arrange(DS)

  SampleOrder <- ds.vec3$Sample
  DR.mat3 <- DR.mat[, match(SampleOrder, colnames(DR.mat))]
  stopifnot(identical(SampleOrder, colnames(DR.mat3)))
  sample.df <- ds.vec3

  col_dii = c("Intact"="mediumaquamarine", "Degraded"="#7C3F11")
  # colwheel <- colortools::wheel("#E41A1C", num=2*2)
  colwheel <- colortools::setColors("#E41A1C", num=2*2)

  column_ha = HeatmapAnnotation(
    "DII"=sample.df$DII,
    "DS"=anno_points(
      sample.df$DS, size=unit(1.2, "mm"),
      axis_param=list(
        at=c(min(sample.df$DS), max(sample.df$DS)),
        labels=sprintf("%.4f", c(min(sample.df$DS), max(sample.df$DS)))
      )
    ),
    "PD"=anno_points(
      sample.df$PD,
      size=unit(1.2, "mm"),
      gp=gpar(col=ifelse(sample.df$PD>cutoff, "#7C3F11", "mediumaquamarine")),
      axis_param=list(
        at=c(min(sample.df$PD), max(sample.df$PD)),
        labels=sprintf("%.4f", c(min(sample.df$PD), max(sample.df$PD)))
      )
    ),
    annotation_name_gp=gpar(fontsize=8),
    annotation_name_rot=c(00),
    col=list(DII=col_dii),
    show_legend=FALSE,
    annotation_name_side=c("left"),
    gap=unit(2, "mm")
  )

  gene.df <- gene.df %>%
    arrange(desc(w_norm), desc(merged_kb))

  GeneOrder <- rownames(gene.df)
  DR.mat4 <- DR.mat3[match(GeneOrder, rownames(DR.mat3)), ]
  stopifnot(identical(GeneOrder, rownames(DR.mat4)))

  mt <- t(gene.df[, c("merged_kb", "mean.TPM")])
  colnames(mt) <- gene.df$geneSymbol

  for (i in 1:nrow(mt)) {
    assign(paste0("mt", i), convert_Cont2eCDF(mt[i,], rownames(mt)[i], margin=1))
  }

  ecdf <- matrix(0, ncol(mt), nrow(mt))
  for (i in 1:nrow(mt)) {
    ecdf[, i] <- get(paste0("mt", i))[, 2]
  }
  rownames(ecdf) <- rownames(mt1)
  colnames(ecdf) <- rownames(mt)

  ecdf <- as.data.frame(ecdf)

  stopifnot(identical(GeneOrder, rownames(ecdf)))
  gene.df <- gene.df %>%
    mutate(merged_kb_ecdf=ecdf$merged_kb,
           mean.TPM_ecdf=ecdf$mean.TPM)

  col_fun_prop3 = colorRamp2(c(0, 0.5, 1), c(colwheel[2], "white", colwheel[4]))
  lgd3 = Legend(col_fun=col_fun_prop3, title="Gene length", title_position="topcenter", at=c(0, 0.5, 1), direction="horizontal", legend_height=unit(0.5, "mm"), legend_width=unit(25, "mm"), labels_gp=gpar(fontsize=8),
                labels=c(round(min(gene.df$merged_kb),1), round(quantile(gene.df$merged_kb, 0.5),1), round(max(gene.df$merged_kb),1)))

  col_fun_prop4 = colorRamp2(c(0, 0.5, 1), c(colwheel[2], "white", colwheel[4]))
  lgd4 = Legend(col_fun=col_fun_prop4, title="Mean TPM", title_position="topcenter", at=c(0, 0.5, 1), direction="horizontal", legend_height=unit(0.5, "mm"), legend_width=unit(25, "mm"), labels_gp=gpar(fontsize=8),
                labels=c(round(min(gene.df$mean.TPM),1), round(quantile(gene.df$mean.TPM, 0.5),1), round(max(gene.df$mean.TPM),1)))

  row_ha = rowAnnotation(
    "Weight"=anno_points(
      gene.df$w_norm,
      size=unit(1.2, "mm"),
      axis_param=list(
        at=c(min(gene.df$w_norm), max(gene.df$w_norm)),
        labels=sprintf("%.4f", c(min(gene.df$w_norm), max(gene.df$w_norm)))
      )
    ),
    "Gene length"=gene.df$merged_kb_ecdf,
    "Mean TPM"=gene.df$mean.TPM_ecdf,
    annotation_name_gp=gpar(fontsize=8),
    annotation_name_rot=c(90),
    col=list(
      "Gene length"=col_fun_prop3,
      "Mean TPM"=col_fun_prop4
    ),
    show_legend=FALSE,
    gap=unit(2, "mm")
  )

  stopifnot(identical(GeneOrder, rownames(DR.mat4)))
  stopifnot(identical(SampleOrder, colnames(DR.mat4)))

  summary <- summary(as.numeric(DR.mat4))
  DR_ul <- round(c(summary[2], summary[5])[which.max(c(summary[2], summary[5]))], 2)

  hm.degrate <- Heatmap(
    DR.mat4,
    name="DR",
    col=colorRamp2(c(-DR_ul, 0, DR_ul), c("blue", "white", "red")),
    top_annotation=column_ha,
    right_annotation=row_ha,
    column_split=factor(ds.vec3$DII, levels=unique(ds.vec3$DII)),
    show_row_names=FALSE,
    show_column_names=FALSE,
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    column_gap=unit(5, "mm"),
    use_raster=TRUE,
    show_heatmap_legend=FALSE,
    width=unit(3.8, "in"),
    height=unit(3.8, "in")
  )

  col_fun_prop0 = colorRamp2(c(-DR_ul, 0, DR_ul), c("blue", "white", "red"))
  lgd0 = Legend(col_fun=col_fun_prop0, title="DR", at=c(-DR_ul, 0, DR_ul), direction="horizontal", legend_height=unit(0.5, "mm"), legend_width=unit(25, "mm"), labels_gp=gpar(fontsize=8),
                labels=c(-DR_ul, 0, DR_ul), title_position="topcenter")

  p <- grid.grabExpr({
    pushViewport(viewport(clip="off")) # to show legends outside of canvas

    draw(hm.degrate)
    pd = packLegend(list=lapply(c(0,3:4), function(i) get(paste0("lgd", i))),
                    max_width=unit(15, "cm"),
                    direction="horizontal", column_gap=unit(1.5, "cm"))
    draw(pd, x=unit(0.45, "npc"), y=unit(0.05, "npc"), just=c("center", "bottom"))
  })

  # Save to outFile
  if (is.null(outFile)) {
    outDir <- file.path(getwd(), "plot")
    if (!dir.exists(outDir)) {
      dir.create(outDir, recursive=TRUE, showWarnings=FALSE)
    }
    outFile <- file.path(outDir, "Fig_DIIwt.png")
  } else {
    outDir <- dirname(outFile)
    if (!dir.exists(outDir)) {
      dir.create(outDir, recursive=TRUE, showWarnings=FALSE)
    }
  }

  message("Saving DIIwt plot to: ", normalizePath(outFile, mustWork=FALSE))

  grDevices::png(file=outFile, width=800*2, height=800, res=100)
  on.exit(dev.off(), add=TRUE)

  dp <- ggpubr::ggarrange(d, as_ggplot(p), nrow=1, ncol=2, align="none")
  print(dp)
}


#' Plot suboptimal/optimal index outputs
#'
#' @param SOIresult outputs from \code{get_SOI} function
#' @param cutoff numeric threshold on projection depth used to classify samples.
#' @param outFile a directory with a file name to save outputs. Default is NULL.
#' @return figures for the distribution of SOI by PD; and the relation of wCV and MCD.
#' @references https://jtr13.github.io/cc21fall2/raincloud-plot-101-density-plot-or-boxplotwhy-not-do-both.html
#' @importFrom dplyr case_when
#' @importFrom ggplot2 ggplot aes position_nudge position_jitter
#' @importFrom ggplot2 geom_point geom_boxplot scale_fill_manual
#' @importFrom ggplot2 scale_color_identity labs guides guide_legend
#' @importFrom ggplot2 theme element_rect element_text coord_flip
#' @importFrom ggplot2 geom_line scale_color_manual xlab ylab coord_cartesian
#' @importFrom ggplot2 margin
#' @export

plot_SOI = function(SOIresult, cutoff=3, outFile=NULL) {

  auc.vec2 <- SOIresult$auc.vec[, c("AUC", "PD")]
  auc.coord <- SOIresult$auc.coord
  rangeMCD <- SOIresult$rangeMCD
  df <- data.frame(var=rep(colnames(auc.vec2), each=nrow(auc.vec2)),
                   value=as.vector(as.matrix(auc.vec2)))
  df$var <- factor(df$var, levels=rev(colnames(auc.vec2)), ordered=TRUE)

  # Distribution of PD by SOI
  make_vioplot <- function(df_sub, cutoff) {
    ggplot2::ggplot(df_sub, aes(x=var, y=value, fill=var)) +
      PupillometryR::geom_flat_violin(position=position_nudge(x=0.2), alpha=0.4) +
      geom_point(aes(color=case_when(
        var=="AUC" ~ "grey",
        var=="PD" & value>cutoff ~ "red",
        var=="PD" & value<=cutoff ~ "darkgreen")),
        position=position_jitter(width=0.15, seed=12345), size=1.5, alpha=0.5, show.legend=F) +
      geom_boxplot(width=0.25, outlier.shape=NA, alpha=0.5) +
      scale_fill_manual(values=c(rep("grey", 2))) +
      scale_color_identity() +
      labs(x="", y="", fill="", title="") +
      guides(fill=guide_legend(nrow=1, byrow=TRUE))+
      theme(legend.position="none",
            legend.margin=margin(-5, 0, 0, 0),
            panel.background=element_rect(fill="gray97"),
            axis.text=element_text(size=25),
            axis.title=ggplot2::element_blank()) +
      coord_flip()
  }

  p1 <- make_vioplot(filter(df, var=="AUC"), cutoff)
  p2 <- make_vioplot(filter(df, var=="PD"), cutoff)

  d <- patchwork::wrap_plots(p1, p2, ncol=1)

  # Relation of wCV and MCD
  auc.coord$SOI <- factor(auc.coord$SOI, levels=c("Suboptimal", "Optimal"))

  p <- ggplot2::ggplot(auc.coord, aes(x=x, y=y, group=Sample, color=SOI)) +
    geom_line(size=1) +
    scale_color_manual(values=c("Suboptimal"=scales::alpha("red", 0.5), "Optimal"=scales::alpha("darkgreen", 0.5))) + # SOI per sample
    geom_rect(data=tibble::tibble(x1=rangeMCD[1], x2=rangeMCD[2], y1=-Inf, y2=+Inf),
              inherit.aes=FALSE,
              mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),
              color="transparent",
              fill="blue",
              alpha=0.07) +
    xlab(expression(paste(log[10], "(MCD+1)"))) + ylab("wCV") +
    coord_cartesian(xlim = c(min(auc.coord$x), (2*rangeMCD[2]-rangeMCD[1]))) +
    guides(size="none") +
    theme(axis.text=element_text(size=25),
          axis.title=element_text(size=25),
          legend.title=element_text(size=18, face="bold"),
          legend.position=c(0.80, 0.85),
          legend.background=element_rect(colour=NA, fill=NA),
          legend.text=element_text(size=18),
          panel.background=element_rect(fill="gray97"))

  # Save to outFile
  if (is.null(outFile)) {
    outDir <- file.path(getwd(), "plot")
    if (!dir.exists(outDir)) {
      dir.create(outDir, recursive=TRUE, showWarnings=FALSE)
    }
    outFile <- file.path(outDir, "Fig_SOI.png")
  } else {
    outDir <- dirname(outFile)
    if (!dir.exists(outDir)) {
      dir.create(outDir, recursive=TRUE, showWarnings=FALSE)
    }
  }

  message("Saving SOI plot to: ", normalizePath(outFile, mustWork=FALSE))

  grDevices::png(file=outFile, width=800*2, height=800, res=100)
  on.exit(dev.off(), add=TRUE)

  dp <- ggpubr::ggarrange(d, p, nrow=1, ncol=2, align="none")
  print(dp)
}
