#' Core helper to compute decay rate
#'
#' @param pileupData exon-only coverage pileup matrix for a single gene.
#' @param exonRanges GRanges object specifying exon coordinates for the gene.
#' @param sampleInfo a sample information table including sample id. The number of rows is equal to the number of samples.
#' @param cases optional character vector specifying a subset of samples.
#'   used for handling missing coverage.
#' @param logshiftVal numeric; passed to \code{process_pileup()}.
#' @param plotNormalization logical; passed to \code{process_pileup()}.
#' @return a numeric vector of decay rates, one value per sample.
#' @details
#' The arguments \code{pileupData}, \code{exonRanges}, \code{logshiftVal}, and
#' \code{plotNormalization} are passed directly to
#' \code{process_pileup()}; see its documentation for details.
#' @references Choi, H.Y., Jo, H., Zhao, X. et al. SCISSOR: a framework for identifying structural changes in RNA transcripts. Nat Commun 12, 286 (2021).
#' @examples
#' ## API illustration only
#' ## Exon-only pileup matrix (rows: exon positions, columns: samples)
#' ## Typically obtained via get_pileupExon()
#' pileupData <- matrix(
#'   c(10, 12, 8,  9,
#'     5,  6,  4,  5),
#'   nrow=2,
#'   byrow=TRUE
#' )
#' colnames(pileupData) <- c("S1", "S2", "S3", "S4")
#'
#' sampleInfo <- data.frame(SampleID=colnames(pileupData))
#'
#' exonRanges <- list(
#'   Gene = "KEAP1",
#'   cRanges = data.frame(
#'     e.start = c(1),
#'     e.end   = c(50001),
#'     row.names = "exon1"
#'   ),
#'   regions     = "chr19:10600000-10650000:+",
#'   new.regions = "chr19:10600000-10650000:+",
#'   strand      = "+"
#' )
#'
#' compute_DR(
#'   pileupData = pileupData,
#'   exonRanges = exonRanges,
#'   sampleInfo = sampleInfo
#' )
#' @export

compute_DR <- function(pileupData, exonRanges, sampleInfo, cases=NULL, logshiftVal=10, plotNormalization=FALSE) {

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
      process_pileup(
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
    allSampleDegRate <- decay.rate.hy(Data=data.process$normalizedData)$slope*d
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
#' @examples
#' ## NOTE:
#' ## This example demonstrates the function interface only.
#' ## Meaningful results require coverage pileup files generated
#' ## from BAM files (see vignette for a full workflow).
#' data("TOY_mrna_mat")
#'
#' ## Interface-only example (no meaningful output is produced)
#' try(
#'   get_DR(
#'     genelist   = TOY_mrna_mat$genes,
#'     pileupPath = rep(NA, length(TOY_mrna_mat$genes)),
#'     sampleInfo = data.frame(SampleID=TOY_mrna_mat$samples),
#'     nCores = 1
#'   ),
#'   silent = TRUE
#' )
#' @export

get_DR <- function(genelist, pileupPath, sampleInfo, cases=NULL, nCores=32) {

  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add=TRUE)

  DR <- foreach::foreach(
    g         = seq_along(pileupPath),
    .combine  = rbind,
    .packages = c("SCISSOR", "RNAshapeQC")
  ) %dopar% {

    # Exon-only coverage for the g-th gene
    pileupData <- get_pileupExon(g, pileupPath, cases)

    if (!is.null(pileupData)&&nrow(pileupData) > 0L) {

      # Exon ranges for this gene
      Gene   <- genelist[g]
      Ranges <- extract_RData(pileupPath[g], "Ranges")
      exonRanges <- get_Ranges(
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
#' @return Invisibly returns NULL; results are saved to \code{outFile}.
#' @references Choi, H.Y., Jo, H., Zhao, X. et al. SCISSOR: a framework for identifying structural changes in RNA transcripts. Nat Commun 12, 286 (2021).
#' @examples
#' ## API illustration only
#' invisible(NULL)
#' @export

gen_DR <- function(Gene, pileupPath, sampleInfo, cases=NULL, Study=NULL, outFile) {

  # Exon-only coverage from multi-study pileupList
  pileupData <- extract_pileupExon(Gene, pileupPath, cases, Study)

  # Gene length based on geneRanges
  geneRanges <- extract_RData(pileupPath, "geneRanges")
  exonRanges <- get_Ranges(
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
#' @examples
#' data("TOY_mrna_mat")
#' get_DIIhc(DR = TOY_mrna_mat$DR)
#' @export

get_DIIhc <- function(DR, topPct=5) {

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


#' Core helper to compute a degraded/intact index using gene weight
#'
#' @param DR a the number of genes x the number of samples matrix of decay rates
#' @param alpha a positive numeric exponent factor to weight the magnitude of decay rates. Default is 2.
#' @param cutoff numeric threshold on projection depth used to classify samples.
#' @param TPM a numeric matrix of TPM values with the same genes in rows and the same samples in columns as \code{DR}.
#' @param thr threshold. Default is 5.
#' @param pct percent. Default is 40.
#' @param genelength a gene length (bp) vector with names as gene IDs.
#' @return a matrix of with decay rate with filtered genes; a matrix including a vector of DII; a data frame of gene info; and a scale factor.
#' @importFrom magrittr %>%
#' @importFrom stats na.omit median uniroot lm coef lowess pnorm
#' @importFrom dplyr mutate select rename inner_join
#' @importFrom graphics abline legend lines text
#' @importFrom grDevices rainbow
#' @importFrom utils data
#' @examples
#' data("TOY_mrna_mat")
#'
#' compute_DIIwt(
#'   DR         = TOY_mrna_mat$DR,
#'   TPM        = TOY_mrna_mat$TPM,
#'   genelength = TOY_mrna_mat$gene_length
#' )
#' @export

compute_DIIwt <- function(DR, alpha=2, cutoff=3, TPM, thr=5, pct=40, genelength) {
  
  if (!identical(rownames(DR), rownames(TPM)))
    stop("DR and TPM must have identical rownames.")
  
  if (!identical(colnames(DR), colnames(TPM)))
    stop("DR and TPM must have identical colnames.")
  
  if (!is.numeric(genelength))
    stop("genelength must be a numeric vector.")
  
  if (is.null(names(genelength)))
    stop("genelength must have gene names.")
  
  if (!identical(names(genelength), rownames(DR)))
    stop("genelength must have same names and order as DR.")
  
  # Remove genes if NAs in DR
  na_rows <- rowSums(is.na(DR)) > 0
  if (any(na_rows)) {
    DR  <- DR[!na_rows, , drop=FALSE]
    TPM <- TPM[!na_rows, , drop=FALSE]
    genelength <- genelength[!na_rows]
  }
  
  # Filter low expressed genes
  genelist2 <- filter_lowExpGenes(genelist=rownames(DR), TPM, thr, pct)

  geneLength <- data.frame(geneid=names(genelength), genelength=genelength) %>%
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
    mutate(PD=pd.rate.hy(DS, qrsc=TRUE), # projection depth
           DII=ifelse(PD>cutoff, "Degraded", "Intact")) # outlier detection

  return(list(DR2=DR2, ds.vec=ds.vec, gene.df=gene.df, s=s))
}


#' Get a degraded/intact index for samples using gene weight
#'
#' @param DR a the number of genes x the number of samples matrix of decay rates
#' @param alpha a positive numeric exponent factor to weight the magnitude of decay rates. Default is 2.
#' @param cutoff numeric threshold on projection depth used to classify samples.
#' @param TPM a numeric matrix of TPM values with the same genes in rows and the same samples in columns as \code{DR}.
#' @param thr threshold. Default is 5.
#' @param pct percent. Default is 40.
#' @param genelength a gene length (bp) vector with names as gene IDs.
#' @param assay.DR character string specifying the assay name containing the DR matrix in a SummarizedExperiment object.
#' @param assay.TPM character string specifying the assay name containing the TPM matrix in a SummarizedExperiment object.
#' @return a matrix of with decay rate with filtered genes; a matrix including a vector of DII; a data frame of gene info; and a scale factor.
#' @importFrom magrittr %>%
#' @importFrom stats na.omit median uniroot
#' @importFrom dplyr mutate select rename inner_join
#' @examples
#' data("TOY_mrna_se")
#'
#' get_DIIwt(TOY_mrna_se)
#' @export

get_DIIwt <- function(DR, alpha=2, cutoff=3, TPM=NULL, thr=5, pct=40, genelength=NULL, assay.DR="DR", assay.TPM="TPM") {
  
  if (inherits(DR, "SummarizedExperiment")) {
    
    se <- DR
    DR  <- SummarizedExperiment::assay(se, assay.DR)
    TPM <- SummarizedExperiment::assay(se, assay.TPM)
    
    if (is.null(genelength)) {
      
      rd <- SummarizedExperiment::rowData(se)
      
      if (!"gene_length" %in% colnames(rd))
        stop("gene_length not found in rowData(se)")
      
      genelength <- rd$gene_length
      names(genelength) <- rownames(se)
    }
  }
  
  if (is.null(TPM))
    stop("TPM must be provided.")
  
  if (is.null(genelength))
    stop("genelength must be provided.")
  
  return(
    compute_DIIwt(
      DR         = DR,
      alpha      = alpha,
      cutoff     = cutoff,
      TPM        = TPM,
      thr        = thr,
      pct        = pct,
      genelength = genelength
    )
  )
}


#' Core helper to compute a mean coverage depth
#'
#' @param pileupData exon-only coverage pileup matrix for a single gene.
#' @param sampleInfo a sample information table including sample id. The number of rows is equal to the number of samples.
#' @param cases optional character vector specifying a subset of samples.
#'   used for handling missing coverage.
#' @return a numeric vector of mean coverage depth, one value per sample.
#' @examples
#' ## API illustration only
#' ## Exon-only pileup matrix (rows: exon positions, columns: samples)
#' ## Typically obtained via get_pileupExon()
#' pileupData <- matrix(
#'   c(10, 12, 8,  9,
#'     5,  6,  4,  5),
#'   nrow=2,
#'   byrow=TRUE
#' )
#' colnames(pileupData) <- c("S1", "S2", "S3", "S4")
#'
#' sampleInfo <- data.frame(SampleID=colnames(pileupData))
#'
#' compute_MCD(
#'   pileupData = pileupData,
#'   sampleInfo = sampleInfo
#' )
#' @export

compute_MCD <- function(pileupData, sampleInfo, cases=NULL) {

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
#' @examples
#' ## NOTE:
#' ## This example demonstrates the function interface only.
#' ## Meaningful results require coverage pileup files generated
#' ## from BAM files (see vignette for a full workflow).
#' data("TOY_total_mat")
#'
#' ## Interface-only example (no meaningful output is produced)
#' try(
#'   get_MCD(
#'     genelist   = TOY_total_mat$genes,
#'     pileupPath = rep(NA, length(TOY_total_mat$genes)),
#'     sampleInfo = data.frame(SampleID=TOY_total_mat$samples),
#'     nCores = 1
#'   ),
#'   silent = TRUE
#' )
#' @export

get_MCD <- function(genelist, pileupPath, sampleInfo, cases=NULL, nCores=32) {

  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add=TRUE)

  MCD <- foreach::foreach(
    g         = seq_along(pileupPath),
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
#' @return Invisibly returns NULL; results are saved to \code{outFile}.
#' @examples
#' ## API illustration only
#' invisible(NULL)
#' @export

gen_MCD <- function(Gene, pileupPath, sampleInfo, cases=NULL, Study=NULL, outFile) {

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
#' @return a numeric vector of window coefficients of variation, one value per sample.
#' @importFrom zoo rollmean rollapply
#' @examples
#' ## API illustration only
#' ## Exon-only pileup matrix (rows: exon positions, columns: samples)
#' ## Typically obtained via get_pileupExon()
#' pileupData <- matrix(
#'   c(10, 12, 8,  9,
#'     5,  6,  4,  5),
#'   nrow=2,
#'   byrow=TRUE
#' )
#' colnames(pileupData) <- c("S1", "S2", "S3", "S4")
#'
#' sampleInfo <- data.frame(SampleID=colnames(pileupData))
#'
#' compute_wCV(
#'   pileupData = pileupData,
#'   sampleInfo = sampleInfo,
#'   rnum       = 10,
#'   winSize    = 2
#' )
#' @export

compute_wCV <- function(pileupData, sampleInfo, rnum=100, method=1, winSize=20, egPct=10, cases=NULL) {

  if (Sys.getenv("_R_CHECK_LIMIT_CORES_", "") != "") {
    return(rep(NA_real_, nrow(sampleInfo)))
  }

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
#' @examples
#' ## NOTE:
#' ## This example demonstrates the function interface only.
#' ## Meaningful results require coverage pileup files generated
#' ## from BAM files (see vignette for a full workflow).
#' data("TOY_total_mat")
#'
#' ## Interface-only example (no meaningful output is produced)
#' try(
#'   get_wCV(
#'     genelist   = TOY_total_mat$genes,
#'     pileupPath = rep(NA, length(TOY_total_mat$genes)),
#'     sampleInfo = data.frame(SampleID=TOY_total_mat$samples),
#'     nCores = 1
#'   ),
#'   silent = TRUE
#' )
#' @export

get_wCV <- function(genelist, pileupPath, sampleInfo, rnum=100, method=1, winSize=20, egPct=10, cases=NULL, nCores=32) {

  if (!(2<=winSize && winSize<=rnum)) {
    stop("The window size ", winSize, " should be in [2, ", rnum, "].")
  }

  cl <- parallel::makeCluster(nCores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add=TRUE)

  wCV <- foreach::foreach(
    g         = seq_along(pileupPath),
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
#' @return Invisibly returns NULL; results are saved to \code{outFile}.
#' @examples
#' ## API illustration only
#' invisible(NULL)
#' @export

gen_wCV <- function(Gene, pileupPath, sampleInfo, rnum=100, method=1, winSize=20, egPct=10, cases=NULL, Study=NULL, outFile) {

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


#' Core helper to compute a suboptimal/optimal index
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
#' @examples
#' data("TOY_total_mat")
#'
#' compute_SOI(
#'   MCD = TOY_total_mat$MCD,
#'   wCV = TOY_total_mat$wCV
#' )
#' @export

compute_SOI <- function(MCD, wCV, rstPct=20, obsPct=50, cutoff=3) {
  
  if (!identical(rownames(MCD), rownames(wCV)))
    stop("MCD and wCV must have identical rownames.")
  
  if (!identical(colnames(MCD), colnames(wCV)))
    stop("MCD and wCV must have identical colnames.")
  
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
    mutate(PD=pd.rate.hy(AUC, qrsc=TRUE), # projection depth
           SOI=ifelse(PD>cutoff, "Suboptimal", "Optimal")) # outlier detection
  
  auc.coord <- smoothData %>%
    select(x, y, Sample) %>%
    inner_join(auc.vec, by="Sample")
  
  return(list(auc.vec=auc.vec, auc.coord=auc.coord, rangeMCD=c(rangeMin, rangeMax)))
}


#' Get a suboptimal/optimal index for samples
#'
#' @param MCD a mean coverage depth is a the number of genes x the number of samples matrix.
#' @param wCV a window coefficient of variation is a the number of genes x the number of samples matrix.
#' @param rstPct restricted percent (one-side) to restrict genes by log transformed MC. Default is 20.
#' @param obsPct span includes the percent of observations in each local regression. Default is 50.
#' @param cutoff numeric threshold on projection depth used to classify samples.
#' @param assay.MCD character string specifying the assay name containing the MCD matrix in a SummarizedExperiment object.
#' @param assay.wCV character string specifying the assay name containing the wCV matrix in a SummarizedExperiment object.
#' @return a matrix including a vector of SOI; a coordinate matrix of smoothed data; and a range of MCD.
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate arrange distinct filter group_by summarise select inner_join
#' @importFrom ggplot2 ggplot_build geom_smooth geom_rect
#' @examples
#' data("TOY_total_se")
#'
#' get_SOI(TOY_total_se)
#' @export

get_SOI <- function(MCD, wCV=NULL, rstPct=20, obsPct=50, cutoff=3, assay.MCD="MCD", assay.wCV="wCV") {

  if (inherits(MCD, "SummarizedExperiment")) {
    
    se <- MCD
    MCD  <- SummarizedExperiment::assay(se, assay.MCD)
    wCV <- SummarizedExperiment::assay(se, assay.wCV)
  }
  
  if (is.null(wCV))
    stop("wCV must be provided.")
  
  return(
    compute_SOI(
      MCD    = MCD,
      wCV    = wCV,
      rstPct = rstPct,
      obsPct = obsPct,
      cutoff = cutoff
    )
  )
}
