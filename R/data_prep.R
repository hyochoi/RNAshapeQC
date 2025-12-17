#' Gene length normalization for a single pileup (for gene 1, sample 1)
#'
#' @param pileup a coverage pileup vector for one sample
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @return the normalized read depth is a vector with length=rnum.
#' @references Choi, H.Y., Jo, H., Zhao, X. et al. SCISSOR: a framework for identifying structural changes in RNA transcripts. Nat Commun 12, 286 (2021).
#' @noRd

norm_pileup.spl <- function(pileup, rnum=100, method=1) {

  row <- c(1:length(pileup))
  depthmat <- cbind(row, pileup)
  pos2 <- data.frame(round(seq(from=1, to=length(pileup), length.out=(2*rnum+1))))
  row_odd <- seq_len(nrow(pos2)) %% 2
  pos3 <- pos2[row_odd==0, ] # even points
  pos4 <- pos2[row_odd==1, ] # odd points
  region <- rep(1:rnum)

  if (method==1) {
    # Method 1: Raw value
    readdepth <- depthmat[pos3,2] # find read depth at even points (green points)

  } else if (method==2) {
    # Method 2: Interpolation
    readdepth <- 10^(round(zoo::rollmean(log10(depthmat[pos4,2]+1),2), 10))-1 # find geometric mean (blue points) using read depth at odd points (red points)

  } else {
    stop(method," is not an option for method.")
  }

  normmat <- as.matrix(cbind(region, pos3, readdepth))
  rownames(normmat) <- region
  colnames(normmat) <- c("region", "pos", "readdepth")

  return(normmat[,c("readdepth")])
}


#' Gene length normalization for a pileup matrix (for gene 1, all samples)
#'
#' @param pileupData a coverage pileup matrix that columns are samples
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @param nCores number of cores used internally for normalization.
#' @return the normalized read depth is a rnum x the number of samples matrix.
#' @references Choi, H.Y., Jo, H., Zhao, X. et al. SCISSOR: a framework for identifying structural changes in RNA transcripts. Nat Commun 12, 286 (2021).
#' @noRd

norm_pileup.gene <- function(pileupData, rnum=100, method=1, nCores=32) {

  if (!(method %in% c(1, 2))) {
    stop(method, " is not an option for method.")
  }

  normmat.gene <- do.call(cbind, parallel::mclapply(seq_len(ncol(pileupData)), function(i) {
    norm_pileup.spl(pileup=pileupData[, i], rnum=rnum, method=method)
  }, mc.cores=max(1L, as.integer(floor(nCores/2)))))
  colnames(normmat.gene) <- colnames(pileupData)

  return(normmat.gene)
}


#' Gene length normalization for pileup lists (for all genes, all samples)
#'
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param geneNames gene names per file. If NULL, Gene i with the same length of pileupPath be set. Default is NULL.
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @param cases a vector of specific samples among all samples in pileup. If NULL, all samples are selected. Default is NULL.
#' @return the normalized read depth is a rnum x the number of samples matrix at each gene list.
#' @references Choi, H.Y., Jo, H., Zhao, X. et al. SCISSOR: a framework for identifying structural changes in RNA transcripts. Nat Commun 12, 286 (2021).
#' @noRd

norm_pileup.list <- function(pileupPath, geneNames=NULL, rnum=100, method=1, cases=NULL) {

  if (is.null(geneNames)) {
    geneNames = paste0("Gene_", c(1:length(pileupPath)))
  }

  geneNames <- geneNames[file.exists(pileupPath)]
  pileupPath <- pileupPath[file.exists(pileupPath)]

  if(length(pileupPath)!=length(geneNames)) {
    stop("pileupPath must be the same length as geneNames")

  } else if(length(pileupPath)==length(geneNames)) {
    pileupList <- list()
    for (g in 1:length(pileupPath)){
      pileupList[[g]] <- list()
      pileupList[[g]] <- get_pileupExon(g, pileupPath, cases)
    }
    if (method==1 | method==2) {
      # Method 1: Raw value
      # Method 2: Interpolation
      normmat.list = lapply(pileupList, FUN=function(x) norm_pileup.gene(pileupData=x, rnum=rnum, method=method))

      # Add gene names in each list
      names(normmat.list) <- geneNames

    } else {
      stop(method," is not an option for method.")
    }
  }

  return(normmat.list)
}


#' Scale normalized transcript coverage using gene length normalization
#'
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param geneNames gene names per file. If NULL, Gene i with the same length of pileupPath be set. Default is NULL.
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @param scale TRUE/FALSE returns the scaled/unscaled normalized transcript coverage. Default is TRUE.
#' @param cases a vector of specific samples among all samples in pileup. If NULL, all samples are selected. Default is NULL.
#' @return gene lists (rows are regions and columns are samples) for the normalized transcript coverage after gene length normalization
#' @references Choi, H.Y., Jo, H., Zhao, X. et al. SCISSOR: a framework for identifying structural changes in RNA transcripts. Nat Commun 12, 286 (2021).
#' @noRd

scale_pileup.list <- function(pileupPath, geneNames=NULL, rnum=100, method=1, scale=TRUE, cases=NULL) {

  # Gene length normalization
  normlist = norm_pileup.list(pileupPath, geneNames, rnum=rnum, method=method, cases=cases)

  # Log-transformation
  log.normlist = lapply(normlist, FUN=function(x) log10(x+1))

  if (is.na(scale) | scale==TRUE) {
    scale.log.normlist = lapply(log.normlist, FUN=function(x) sweep(x, 2, apply(x,2,sum)+0.01, FUN="/"))
    return(scale.log.normlist)

  } else if (scale==FALSE) {
    return(log.normlist)

  } else {
    stop(scale," is not an option for scale.")
  }
}


#' Get metrics from scaled normalized transcript coverage
#'
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param geneNames gene names per file. If NULL, Gene i with the same length of pileupPath be set. Default is NULL.
#' @param rnum the number of regions for uniformly dividing the x-axis. Default is 100.
#' @param method 1 and 2 return the raw read depth and the interpolated read depth at the normalized genomic position, respectively. Default is 1.
#' @param scale TRUE/FALSE returns the scaled/unscaled normalized transcript coverage. Default is TRUE.
#' @param margin 1, 2, and 3 return metrics per sample, per gene, and across the genes per sample, respectively.
#' @return metrics including mean, sd, CV (sd/mean), median, mad, and robustCV (mad/median) per margin
#' @references Choi, H.Y., Jo, H., Zhao, X. et al. SCISSOR: a framework for identifying structural changes in RNA transcripts. Nat Commun 12, 286 (2021).
#' @importFrom dplyr desc
#' @importFrom stats sd mad quantile
#' @noRd

get_metrics <- function(pileupPath, geneNames=NULL, rnum=100, method=1, scale=TRUE, margin) {

  scale.log.normlist = scale_pileup.list(pileupPath, geneNames, rnum=rnum, method=method, scale=TRUE)
  mar <- list(2, 3, 2:3)
  array_data <- simplify2array(scale.log.normlist)

  if (margin %in% 1:3) {
    # Use positive values only; if all values are 0 then all stats are 0
    pos_mask <- array_data>0
    var.sum <- apply(array_data, mar[[margin]], function(x) sum(x^2))
    pos_data <- replace(array_data, !pos_mask, NA)  # replace non-positive values with NA
    var.mean <- apply(pos_data, mar[[margin]], mean, na.rm=TRUE)
    var.sd <- apply(pos_data, mar[[margin]], sd, na.rm=TRUE)
    var.median <- apply(pos_data, mar[[margin]], median, na.rm=TRUE)
    var.mad <- apply(pos_data, mar[[margin]], mad, na.rm=TRUE)
    CV <- ifelse(var.mean<1e-10 | var.sd<1e-10, 0, var.sd/var.mean) # to adjust NaN, Inf
    robustCV <- ifelse(var.median<1e-10 | var.mad<1e-10, 0, var.mad/var.median)

  } else {
    stop(margin, " is not an option for margin.")
  }

  # Assemble the metrics
  if (margin %in% 1:2) {
    metrics <- data.frame(
      mean=var.mean,
      sd=var.sd,
      CV=CV,
      median=var.median,
      mad=var.mad,
      robustCV=robustCV
    )

  } else if (margin==3) {
    metrics <- data.frame(
      convert_pivot.longer(var.mean, c("sample", "gene", "mean")),
      sd=convert_pivot.longer(var.sd, c("sample", "gene", "sd"))[, "sd"],
      CV=convert_pivot.longer(CV, c("sample", "gene", "CV"))[, "CV"],
      median=convert_pivot.longer(var.median, c("sample", "gene", "median"))[, "median"],
      mad=convert_pivot.longer(var.mad, c("sample", "gene", "mad"))[, "mad"],
      robustCV=convert_pivot.longer(robustCV, c("sample", "gene", "robustCV"))[, "robustCV"]
    )
  }

  return(metrics)
}
