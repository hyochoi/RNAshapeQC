#' Core helper to compute decay rate
#'
#' @param pileupData exon-only coverage pileup matrix for a single gene.
#' @param exonRanges GRanges object specifying exon coordinates for the gene.
#' @param sampleInfo a sample information table including sample id. The number of rows is equal to the number of samples.
#' @param cases optional character vector specifying a subset of samples.
#'   used for handling missing coverage.
#' @param logshiftVal numeric; passed to \code{SCISSOR::process_pileup()}.
#' @param plotNormalization logical; passed to \code{SCISSOR::process_pileup()}.
#' @details
#' The arguments \code{pileupData}, \code{exonRanges}, \code{logshiftVal}, and
#' \code{plotNormalization} are passed directly to
#' \code{SCISSOR::process_pileup()}; see its documentation for details.
#' @references Choi, H.Y., Jo, H., Zhao, X. et al. SCISSOR: a framework for identifying structural changes in RNA transcripts. Nat Commun 12, 286 (2021).
#' @export

compute_DR = function(pileupData, exonRanges, sampleInfo, cases=NULL, logshiftVal=10, plotNormalization=FALSE) {

  # No data: return NA vector
  if (is.null(pileupData)||nrow(pileupData) == 0L) {
    if (is.null(cases)) {
      return(rep(NA_real_, nrow(sampleInfo)))
    } else {
      return(rep(NA_real_, length(cases)))
    }
  }

  d <- nrow(pileupData)

  # Run process_pileup with basic error handling
  data.process <- tryCatch(
    {
      SCISSOR::process_pileup(
        pileupData        = pileupData,
        Ranges            = exonRanges,
        logshiftVal       = logshiftVal,
        plotNormalization = plotNormalization
      )
    },
    error = function(e) {
      message("process_pileup error: ", conditionMessage(e))
      return(NULL)
    }
  )

  # If process_pileup failed, return NA vector
  if (is.null(data.process)) {
    if (is.null(cases)) {
      allSampleDegRate <- rep(NA_real_, nrow(sampleInfo))
    } else {
      allSampleDegRate <- rep(NA_real_, length(cases))
    }
  } else {
    allSampleDegRate <- SCISSOR::decay.rate.hy(Data=data.process$normalizedData)$slope*d
  }

  return(allSampleDegRate)
}


#' Get a decay rate for genes and samples (for a genelist)
#'
#' @param genelist a vector of gene names
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param sampleInfo a sample information table including sample id. The number of rows is equal to the number of samples.
#' @param cases a vector of specific samples among all samples in pileup. If NULL, all samples are selected. Default is NULL.
#' @param nCores the number of cores for parallel computing. Default is 32.
#' @return DR is a the number of genes x the number of samples matrix.
#' @references Choi, H.Y., Jo, H., Zhao, X. et al. SCISSOR: a framework for identifying structural changes in RNA transcripts. Nat Commun 12, 286 (2021).
#' @importFrom foreach %dopar%
#' @export

get_DR = function(genelist, pileupPath, sampleInfo, cases=NULL, nCores=32) {

  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add=TRUE)

  DR <- foreach::foreach(
    g         = 1:length(pileupPath),
    .combine  = rbind,
    .packages = c("SCISSOR", "RNAshapeQC")
  ) %dopar% {

    # Exon-only coverage for the g-th gene
    pileupData <- get_pileupExon(g, pileupPath, cases)

    if (!is.null(pileupData)&&nrow(pileupData) > 0L) {

      # Exon ranges for this gene
      Gene   <- genelist[g]
      Ranges <- extract_RData(pileupPath[g], "Ranges")
      exonRanges <- SCISSOR::get_Ranges(
        Gene       = Gene,
        regions    = Ranges$regions,
        outputType = "only_exon"
      )

      # Compute DR via core helper
      allSampleDegRate <- compute_DR(
        pileupData        = pileupData,
        exonRanges        = exonRanges,
        sampleInfo        = sampleInfo,
        cases             = cases,
        logshiftVal       = 10,
        plotNormalization = FALSE
      )

      # Keep original sample order
      names(allSampleDegRate) <- colnames(pileupData)
      allSampleDegRate

    } else {
      # No coverage for this gene: NA for all samples
      rep(NA_real_, nrow(sampleInfo))
    }
  }

  rownames(DR) <- genelist
  return(DR)
}


#' Get a decay rate for genes and samples (for a single gene)
#'
#' @param Gene a character of gene name
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param sampleInfo a sample information table including sample id. The number of rows is equal to the number of samples.
#' @param cases a vector of specific samples among all samples in pileup. If NULL, all samples are selected. Default is NULL.
#' @param Study a character of study abbreviation in the pileupList. Default is NULL.
#' @param outFile a directory with a file name to save outputs
#' @return DR is a the number of genes x the number of samples matrix.
#' @references Choi, H.Y., Jo, H., Zhao, X. et al. SCISSOR: a framework for identifying structural changes in RNA transcripts. Nat Commun 12, 286 (2021).
#' @export

gen_DR = function(Gene, pileupPath, sampleInfo, cases=NULL, Study=NULL, outFile) {

  # Exon-only coverage from multi-study pileupList
  pileupData <- extract_pileupExon(Gene, pileupPath, cases, Study)

  # Gene length based on geneRanges
  geneRanges <- extract_RData(pileupPath, "geneRanges")
  exonRanges <- SCISSOR::get_Ranges(
    Gene       = Gene,
    regions    = geneRanges$regions,
    outputType = "only_exon"
  )
  genelength <- max(exonRanges$cRanges)

  # Compute DR via core helper
  allSampleDegRate <- compute_DR(
    pileupData        = pileupData,
    exonRanges        = exonRanges,
    sampleInfo        = sampleInfo,
    cases             = cases,
    logshiftVal       = 10,
    plotNormalization = FALSE
  )

  # Save outputs
  save(genelength, allSampleDegRate, file=outFile)
}


#' Get a degraded/intact index for samples using hierarchical clustering
#'
#' @param DR a the number of genes x the number of samples matrix of decay rates
#' @param topPct top percentages of decay rates defined as degrateGrp=1. Default is 5.
#' @return a matrix of binary converted decay rates; hierarchical clustering outputs of samples; and a vector of DII per sample.
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise arrange mutate left_join select
#' @importFrom stats na.omit dist hclust as.dendrogram
#' @export

get_DIIhc = function(DR, topPct=5) {

  DR <- na.omit(DR)
  # degrateGrp <- ifelse(DR>=quantile(as.vector(DR), 1-topPct/100, na.rm=TRUE), 1, 0)
  degrateGrp <- 1*(DR>=stats::quantile(as.vector(DR), 1-topPct/100, na.rm=TRUE))
  hc_columns <- hclust(dist(t(degrateGrp), method="manhattan"), method="ward.D")

  # Convert to dendrogram
  dend_columns <- as.dendrogram(hc_columns)

  hd.vec <- data.frame(Sample=hc_columns$labels,
                       Cluster=dendextend::cutree(dend_columns, k=2),
                       Order=hc_columns$order,
                       ColSums=colSums(degrateGrp))

  # Define index by meanSum
  cluster_means <- hd.vec %>%
    group_by(Cluster) %>%
    summarise(meanSum=base::mean(ColSums), .groups="drop") %>%
    arrange(desc(meanSum)) %>%
    mutate(DII=c("Degraded", "Intact"), # Degraded if meanSum is greater
           Cluster_new=1:2) # Degraded: Cluster=1

  hd.vec <- hd.vec %>%
    left_join(cluster_means, by="Cluster") %>%
    select(Sample, Cluster=Cluster_new, Order, ColSums, DII)

  return(list(degrateGrp=degrateGrp, hc_columns=hc_columns, hd.vec=hd.vec, cluster_means=cluster_means))
}


#' Get a degraded/intact index for samples using gene weight
#'
#' @param DR a the number of genes x the number of samples matrix of decay rates
#' @param alpha a positive numeric exponent factor to weight the magnitude of decay rates. Default is 2.
#' @param cutoff numeric threshold on projection depth used to classify samples.
#' @param TPM a numeric matrix of TPM values with the same genes in rows and the same samples in columns as \code{DR}.
#' @param thru threshold. Default is 5.
#' @param pct percent. Default is 40.
#' @param genelength.mat a one-column gene length (bp) matrix with row names as gene IDs (no column names).
#' @return a matrix of with decay rate with filtered genes; a matrix including a vector of DII; a data frame of gene info; and a scale factor.
#' @importFrom magrittr %>%
#' @importFrom stats na.omit median uniroot
#' @importFrom dplyr mutate select rename inner_join
#' @export

get_DIIwt = function(DR, alpha=2, cutoff=3, TPM, thru=5, pct=40, genelength.mat) {

  if (!identical(colnames(DR), colnames(TPM))) {
    stop("DR and TPM should have the same samples.")
  }

  DR <- na.omit(DR)
  genelist2 <- filter_lowExpGenes(genelist=rownames(DR), TPM, thru, pct)

  geneLength <- data.frame(geneid=rownames(genelength.mat), genelength=genelength.mat) %>%
    mutate(merged_kb=genelength/1000) %>%
    select(geneid, merged_kb)

  mean.TPM <- apply(TPM, 1, function(x) mean(x, na.rm=TRUE))
  TPMdf <- data.frame(Gene=names(mean.TPM), mean.TPM=mean.TPM)

  # Scale factor
  gene.df <- geneLength[match(genelist2, geneLength$geneid), ] %>%
    rename("geneSymbol"="geneid") %>%
    inner_join(TPMdf, by=c("geneSymbol"="Gene")) %>%
    mutate(logTPM=log2(mean.TPM+1))

  f <- function(s) {
    median(log2(gene.df$merged_kb*1000/s+1)) - median(gene.df$logTPM)
  }

  sol <- uniroot(f, interval=c(1e-6, 1e6))
  s <- sol$root

  # Gene weight using mean TPM and gene length
  gene.df <- gene.df %>%
    mutate(scaledLength=merged_kb*1000/s,
           logLength=log2(scaledLength+1),
           w_raw=logTPM*logLength,
           w_norm=w_raw/sum(w_raw)) %>% # weight normalized to sum 1
    select(geneSymbol, merged_kb, mean.TPM, scaledLength, logTPM, logLength, w_norm)

  rownames(gene.df) <- gene.df$geneSymbol

  # Update with filtered genes
  DR2 <- DR[match(genelist2, rownames(DR)), ]

  # Review
  stopifnot(identical(rownames(DR2), gene.df$geneSymbol)) # gene order
  stopifnot(abs(sum(gene.df$w_norm)-1) < 1e-8) # sum 1

  # Signed degradation score
  w_norm=gene.df$w_norm
  ds.vec <- data.frame(Sample=colnames(DR2),
                       DS=as.vector(crossprod(w_norm, sign(DR2)*abs(DR2)^alpha))) %>%
    mutate(PD=SCISSOR::pd.rate.hy(DS, qrsc=TRUE), # projection depth
           DII=ifelse(PD>cutoff, "Degraded", "Intact")) # outlier detection

  return(list(DR2=DR2, ds.vec=ds.vec, gene.df=gene.df, s=s))
}


#' Core helper to compute a mean coverage depth
#'
#' @param pileupData exon-only coverage pileup matrix for a single gene.
#' @param sampleInfo a sample information table including sample id. The number of rows is equal to the number of samples.
#' @param cases optional character vector specifying a subset of samples.
#'   used for handling missing coverage.
#' @export

compute_MCD = function(pileupData, sampleInfo, cases=NULL) {

  # No data: return NA vector
  if (is.null(pileupData)||nrow(pileupData) == 0L) {
    if (is.null(cases)) {
      return(rep(NA_real_, nrow(sampleInfo)))
    } else {
      return(rep(NA_real_, length(cases)))
    }
  }

  # Use positive values only; if all values are 0 then all stats are 0
  sum_cov <- colSums(pmax(pileupData, 0), na.rm=TRUE)
  n_pos   <- colSums(pileupData>0, na.rm=TRUE)
  mcd <- ifelse(sum_cov<1e-10|n_pos<1e-10, 0, sum_cov/n_pos) # to adjust NaN, Inf

  return(mcd)
}


#' Get a mean coverage depth for genes and samples (for a genelist)
#'
#' @param genelist a vector of gene names
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param sampleInfo a sample information table including sample id. The number of rows is equal to the number of samples.
#' @param cases a vector of specific samples among all samples in pileup. If NULL, all samples are selected. Default is NULL.
#' @param nCores the number of cores for parallel computing. Default is 32.
#' @return MCD is a the number of genes x the number of samples matrix.
#' @importFrom foreach %dopar%
#' @export

get_MCD = function(genelist, pileupPath, sampleInfo, cases=NULL, nCores=32) {

  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add=TRUE)

  MCD <- foreach::foreach(
    g         = 1:length(pileupPath),
    .combine  = rbind,
    .packages = c("SCISSOR", "RNAshapeQC")
  ) %dopar% {

    pileupData <- get_pileupExon(g, pileupPath, cases)

    if (!is.null(pileupData)&&nrow(pileupData) > 0L) {

      mcd.vec <- compute_MCD(
        pileupData = pileupData,
        sampleInfo = sampleInfo,
        cases      = cases
      )

      # Keep original sample ordering / names when available
      if (!is.null(colnames(pileupData))) {
        names(mcd.vec) <- colnames(pileupData)
      }
      mcd.vec

    } else {
      # No coverage for this gene: NA for all samples
      rep(NA_real_, nrow(sampleInfo))
    }
  }

  rownames(MCD) <- genelist
  return(MCD)
}


#' Get a mean coverage depth for genes and samples (for a single gene)
#'
#' @param Gene a character of gene name
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param sampleInfo a sample information table including sample id. The number of rows is equal to the number of samples.
#' @param cases a vector of specific samples among all samples in pileup. If NULL, all samples are selected. Default is NULL.
#' @param Study a character of study abbreviation in the pileupList. Default is NULL.
#' @param outFile a directory with a file name to save outputs
#' @return MCD is a the number of genes x the number of samples matrix.
#' @export

gen_MCD = function(Gene, pileupPath, sampleInfo, cases=NULL, Study=NULL, outFile) {

  # Exon-only coverage from single-study pileup or multi-study pileupList
  pileupData <- extract_pileupExon(Gene, pileupPath, cases, Study)

  # Compute MCD via core helper
  allSampleMCD <- compute_MCD(
    pileupData = pileupData,
    sampleInfo = sampleInfo,
    cases      = cases
  )

  save(allSampleMCD, file=outFile)
}


#' Core helper to compute a window coefficient of variation
#'
#' @param pileupData exon-only coverage pileup matrix for a single gene.
#' @param sampleInfo a sample information table including sample id. The number of rows is equal to the number of samples.
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @param winSize window size of the rolling window. Default is 20.
#' @param egPct edge percent (one-side) to calculate the trimmed mean. Default is 10.
#' @param cases optional character vector specifying a subset of samples.
#'   used for handling missing coverage.
#' @importFrom zoo rollmean rollapply
#' @export

compute_wCV = function(pileupData, sampleInfo, rnum=100, method=1, winSize=20, egPct=10, cases=NULL) {

  # No data: return NA vector
  if (is.null(pileupData)||nrow(pileupData) == 0L) {
    if (is.null(cases)) {
      return(rep(NA_real_, nrow(sampleInfo)))
    } else {
      return(rep(NA_real_, length(cases)))
    }
  }

  # Gene-level normalization
  norm_pileup <- norm_pileup.gene(pileupData, rnum=rnum, method=method)

  # Rolling CV for each sample
  rmean <- zoo::rollmean(norm_pileup, winSize, fill=NA, align="center", by.column=TRUE)
  rsd <- zoo::rollapply(norm_pileup, winSize, stats::sd, fill=NA, align="center", by.column=TRUE)
  cv.mat <- ifelse(rmean<1e-10 | rsd<1e-10, 0, rsd/rmean) # to adjust NaN, Inf

  # 0-adjusted trimmed mean per sample
  trmean_col <- function(column, trimFrac) {
    posVals <- column[column>0]
    if (all(is.na(posVals))) {
      return(NA_real_)
    } else {
      return(base::mean(posVals, na.rm=TRUE, trim=trimFrac))
    }
  }
  trimFrac <- egPct/100
  wcv.vec <- apply(cv.mat, 2, trmean_col, trimFrac=trimFrac)

  return(wcv.vec)
}


#' Get a window coefficient of variation for genes and samples (for a genelist)
#'
#' @param genelist a vector of gene names
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param sampleInfo a sample information table including sample id. The number of rows is equal to the number of samples.
#' @param rnum the number of regions for uniformly dividing the x-axis for gene length normalization. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @param winSize window size of the rolling window. Default is 20.
#' @param egPct edge percent (one-side) to calculate the trimmed mean. Default is 10.
#' @param cases a vector of specific samples among all samples in pileup. If NULL, all samples are selected. Default is NULL.
#' @param nCores the number of cores for parallel computing. Default is 32.
#' @return wCV is a the number of genes x the number of samples matrix.
#' @importFrom foreach %dopar%
#' @export

get_wCV = function(genelist, pileupPath, sampleInfo, rnum=100, method=1, winSize=20, egPct=10, cases=NULL, nCores=32) {

  if (!(2<=winSize && winSize<=rnum)) {
    stop("The window size ", winSize, " should be in [2, ", rnum, "].")
  }

  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add=TRUE)

  wCV <- foreach::foreach(
    g         = 1:length(pileupPath),
    .combine  = rbind,
    .packages = c("zoo", "RNAshapeQC")
  ) %dopar% {

    pileupData <- get_pileupExon(g, pileupPath, cases)

    if (!is.null(pileupData)&&nrow(pileupData) > 0L) {

      wcv.vec <- compute_wCV(
        pileupData = pileupData,
        sampleInfo = sampleInfo,
        rnum       = rnum,
        method     = method,
        winSize    = winSize,
        egPct      = egPct,
        cases      = cases
      )

      if (!is.null(colnames(pileupData))) {
        names(wcv.vec) <- colnames(pileupData)
      }
      wcv.vec

    } else {
      # No coverage for this gene: NA for all samples
      rep(NA_real_, nrow(sampleInfo))
    }
  }

  rownames(wCV) <- genelist
  return(wCV)
}


#' Get a window coefficient of variation for genes and samples (for a single gene)
#'
#' @param Gene a character of gene name
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param sampleInfo a sample information table including sample id. The number of rows is equal to the number of samples.
#' @param rnum the number of regions for uniformly dividing the x-axis for gene length normalization. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @param winSize window size of the rolling window. Default is 20.
#' @param egPct edge percent (one-side) to calculate the trimmed mean. Default is 10.
#' @param cases a vector of specific samples among all samples in pileup. If NULL, all samples are selected. Default is NULL.
#' @param Study a character of study abbreviation in the pileupList. Default is NULL.
#' @param outFile a directory with a file name to save outputs
#' @return wCV is a the number of genes x the number of samples matrix.
#' @export

gen_wCV = function(Gene, pileupPath, sampleInfo, rnum=100, method=1, winSize=20, egPct=10, cases=NULL, Study=NULL, outFile) {

  if (!(2<=winSize && winSize<=rnum)) {
    stop("The window size ", winSize, " should be in [2, ", rnum, "].")
  }

  # Exon-only coverage from single-study pileup or multi-study pileupList
  pileupData <- extract_pileupExon(Gene, pileupPath, cases, Study)

  # Compute wCV via core helper
  allSamplewCV <- compute_wCV(
    pileupData = pileupData,
    sampleInfo = sampleInfo,
    rnum       = rnum,
    method     = method,
    winSize    = winSize,
    egPct      = egPct,
    cases      = cases
  )

  save(allSamplewCV, file=outFile)
}


#' Get a suboptimal/optimal index for samples
#'
#' @param MCD a mean coverage depth is a the number of genes x the number of samples matrix.
#' @param wCV a window coefficient of variation is a the number of genes x the number of samples matrix.
#' @param rstPct restricted percent (one-side) to restrict genes by log transformed MC. Default is 20.
#' @param obsPct span includes the percent of observations in each local regression. Default is 50.
#' @param cutoff numeric threshold on projection depth used to classify samples.
#' @return a matrix including a vector of SOI; a coordinate matrix of smoothed data; and a range of MCD.
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate arrange distinct filter group_by summarise select inner_join
#' @importFrom ggplot2 ggplot_build geom_smooth geom_rect
#' @export

get_SOI = function(MCD, wCV, rstPct=20, obsPct=50, cutoff=3) {

  auc.coord <- stats::na.omit(data.frame(Gene=rep(rownames(MCD), ncol(MCD)),
                                  Sample=rep(colnames(MCD), each=nrow(MCD)),
                                  MCD=as.vector(MCD),
                                  wCV=as.vector(wCV))) %>%
    mutate(xMCD=log10(MCD+1)) %>%
    arrange(Sample, xMCD) # sort x-points for AUC

  # LOESS regression
  p <- ggplot2::ggplot(auc.coord, aes(x=xMCD, y=wCV)) +
    geom_smooth(data=auc.coord, aes(group=Sample), method="loess", span=obsPct/100, se=FALSE)
  smoothData <- ggplot_build(p)$data[[1]]

  # Map back to original group
  group_mapping <- auc.coord %>%
    distinct(Sample)
  group_mapping$group_id <- as.numeric(factor(group_mapping$Sample))
  smoothData <- smoothData %>%
    mutate(Sample=group_mapping$Sample[group])

  # Range of MCD
  posMCD <- MCD[MCD>0]
  rangeMin = log10(stats::quantile(posMCD, probs=rstPct/100, na.rm=TRUE)+1)
  rangeMax = log10(stats::quantile(posMCD, probs=1-rstPct/100, na.rm=TRUE)+1)

  auc.vec <- smoothData %>%
    filter(x>=rangeMin & x<rangeMax) %>% # restricted MCD
    group_by(Sample) %>%
    summarise(AUC=DescTools::AUC(x, y, method="spline")) %>% # calculate AUC
    mutate(PD=SCISSOR::pd.rate.hy(AUC, qrsc=TRUE), # projection depth
           SOI=ifelse(PD>cutoff, "Suboptimal", "Optimal")) # outlier detection

  auc.coord <- smoothData %>%
    select(x, y, Sample) %>%
    inner_join(auc.vec, by="Sample")

  return(list(auc.vec=auc.vec, auc.coord=auc.coord, rangeMCD=c(rangeMin, rangeMax)))
}
