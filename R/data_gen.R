#' Internal helper to build regions string for a gene
#'
#' @noRd

.build_regionsStr <- function(Gene, regionsFile, regionsFormat=c("auto", "SCISSOR_gaf", "gencode.regions"), geneCol=1, regionsCol=NULL) {

  regionsFormat <- match.arg(regionsFormat)

  # Read regionsFile (path or data.frame)
  if (is.character(regionsFile)&&length(regionsFile) == 1L) {
    df <- tryCatch(
      utils::read.table(regionsFile, header=TRUE, sep="\t",
                        stringsAsFactors=FALSE, quote="", comment.char=""),
      error = function(e) {
        utils::read.table(regionsFile, header=FALSE, sep="\t",
                          stringsAsFactors=FALSE, quote="", comment.char="")
      }
    )
  } else if (is.data.frame(regionsFile)) {
    df <- regionsFile
  } else {
    stop("`regionsFile` must be a file path or a data.frame.")
  }

  # Auto-detect format if needed
  if (regionsFormat=="auto") {
    if (all(c("gene_name", "regions") %in% colnames(df))) {
      regionsFormat <- "SCISSOR_gaf"
    } else {
      regionsFormat <- "gencode.regions"
    }
  }

  # Extract regions string by format
  if (regionsFormat=="SCISSOR_gaf") {
    if (!all(c("gene_name", "regions") %in% colnames(df))) {
      stop("regionsFormat='SCISSOR_gaf' but 'gene_name' / 'regions' columns not found.")
    }
    if (!Gene %in% df$gene_name) {
      stop(Gene, " is not in `gene_name` of regionsFile.")
    }
    regions <- df$regions[df$gene_name == Gene][1]

  } else if (regionsFormat=="gencode.regions") {
    if (is.null(regionsCol)) regionsCol <- 2L
    if (!Gene %in% df[, geneCol]) {
      stop(Gene, " is not found in column ", geneCol, " of regionsFile.")
    }
    regions <- df[df[, geneCol]==Gene, regionsCol][1]

  } else {
    stop("Unsupported regionsFormat: ", regionsFormat)
  }

  as.character(regions)
}


#' Core helper to construct a pileup from BAM files (for single-study)
#'
#' @noRd

.construct_pileupStudy <- function(Gene, regions, Ranges, BAMfiles, caseIDs, max_depth=100000, strand.specific=FALSE, nCores=10) {

  if (length(BAMfiles)!=length(caseIDs)) {
    stop("`BAMfiles` and `caseIDs` must have the same length.")
  }

  if (is.null(names(BAMfiles))) {
    names(BAMfiles) <- caseIDs
  }

  # Helper: read one BAM file and return coverage vector for the given regions
  read_aBAM <- function(BAM, regions, strand.specific=FALSE, max_depth=100000, ...) {
    bf <- Rsamtools::BamFile(BAM)

    # Parse regions string into GRanges
    chr     <- strsplit(regions, ":")[[1]][1]
    strtend <- do.call(rbind, strsplit(strsplit(strsplit(regions, ":")[[1]][2], ",")[[1]], "-"))
    strnd   <- strsplit(regions, ":")[[1]][3]
    gr      <- GenomicRanges::GRanges(chr, IRanges::IRanges(start=as.numeric(strtend[, 1]), end=as.numeric(strtend[, 2])), strnd)

    # https://www.rdocumentation.org/packages/Rsamtools/versions/1.24.0/topics/pileup
    s_param <- Rsamtools::ScanBamParam(which=gr, what=c("pos"))
    p_param <- Rsamtools::PileupParam(
      distinguish_strands     = strand.specific,
      distinguish_nucleotides = FALSE,
      include_deletions       = FALSE,
      include_insertions      = FALSE,
      left_bins  = NULL,
      query_bins = NULL,
      cycle_bins = NULL,
      max_depth  = max_depth,
      ...
    )

    res <- Rsamtools::pileup(
      bf,
      scanBamParam = s_param,
      pileupParam  = p_param
    )

    # If strand-specific, keep only "+" strand rows if available
    if (strand.specific && "strand" %in% colnames(res)) {
      res <- res[res$strand=="+", , drop=FALSE]
    }

    # Fill zero-coverage positions explicitly
    strtend.num <- matrix(as.numeric(strtend), ncol=2)
    allPos <- unlist(
      lapply(
        seq_len(nrow(strtend.num)),
        function(i) strtend.num[i, 1]:strtend.num[i, 2]
      )
    )

    if (nrow(res)>0L) {
      tmpDepth <- rbind(
        res[, c("pos", "count")],
        cbind(
          pos   = allPos[!allPos %in% res$pos],
          count = rep(0L, length(allPos[!allPos %in% res$pos]))
        )
      )
    } else {
      tmpDepth <- cbind(
        pos   = allPos,
        count = rep(0L, length(allPos))
      )
    }

    outpileup <- tmpDepth[order(tmpDepth[, "pos"]), , drop=FALSE]
    rownames(outpileup) <- seq_len(nrow(outpileup))
    outpileup[, "count"]
  }

  # Use Ranges$new.regions to define the x-axis for the pileup
  new.regions <- Ranges$new.regions
  strtend     <- do.call(rbind, strsplit(strsplit(strsplit(new.regions, ":")[[1]][2], ",")[[1]], "-"))
  strnd       <- strsplit(new.regions, ":")[[1]][3]
  strtend.num <- matrix(as.numeric(strtend), ncol=2)
  allPos      <- unlist(
    lapply(
      seq_len(nrow(strtend.num)),
      function(i) strtend.num[i, 1]:strtend.num[i, 2]
    )
  )

  # Compute coverage for each BAM file in parallel
  pileup <- parallel::mclapply(
    seq_along(BAMfiles),
    FUN = function(i) {
      bam <- BAMfiles[i]
      message("Reading BAM: ", bam)
      read_aBAM(
        BAM            = bam,
        regions        = new.regions,
        strand.specific = strand.specific,
        max_depth      = max_depth
      )
    },
    mc.cores = nCores
  ) |>
    do.call(what = "cbind")

  rownames(pileup) <- allPos
  colnames(pileup) <- caseIDs

  # Reverse rows if the strand is "-"
  if (strnd=="-") {
    pileup <- pileup[rev(seq_len(nrow(pileup))), , drop=FALSE]
  }

  list(pileup=pileup, regions=regions, Ranges=Ranges)
}


#' Construct a per-gene pileup from BAM files (for single-study or multi-study)
#'
#' @param Gene a character of gene name
#' @param studylist a character vector of study IDs or abbreviation/name
#' @param regionsFile either a file path or a data.frame specifying gene regions.
#'   If a file path is given, it is read using \code{read.table()} with automatic
#'   handling of header/non-header cases.
#' @param regionsFormat character; one of \code{"auto"}, \code{"SCISSOR_gaf"},
#'   or \code{"gencode.regions"}. In \code{"auto"} mode, the function attempts
#'   to detect whether the file uses SCISSOR-style columns \code{gene_name} and
#'   \code{regions}, falling back to \code{"gencode.regions"} otherwise.
#' @param geneCol integer; column index for the gene identifier when
#'   \code{regionsFormat = "gencode.regions"}. Default is 1.
#' @param regionsCol integer; column index for the regions string when
#'   \code{regionsFormat = "gencode.regions"}. Default is \code{NULL}, which
#'   is interpreted as 2.
#' @param bamFilesList named list of character vectors of BAM file paths
#' @param caseIDList named list of character vectors of sample IDs corresponding
#'   to \code{bamFilesList}
#' @param max_depth integer; max depth parameter for \code{Rsamtools::PileupParam}.
#'   Default is 100000.
#' @param strand.specific Logical; whether to use strand-specific pileup.
#'   If \code{TRUE}, only the "+" strand is retained (when strand information
#'   is available). Default is \code{FALSE}.
#' @param nCores the number of cores for parallel computing. Default is 10.
#' @param outFile a directory with a file name to save outputs. Default is NULL.
#' @return a pileup matrix, regions, and ranges of genomic positions
#' @examples
#' ## API illustration only
#' invisible(NULL)
#' @export

construct_pileup <- function(
    Gene,
    studylist,
    regionsFile,
    regionsFormat = c("auto", "SCISSOR_gaf", "gencode.regions"),
    geneCol      = 1,
    regionsCol   = NULL,
    bamFilesList,
    caseIDList,
    max_depth       = 100000,
    strand.specific = FALSE,
    nCores          = 10,
    outFile         = NULL
) {

  regionsFormat <- match.arg(regionsFormat)

  # Build regions string for this gene from the annotation file / data.frame
  regions <- .build_regionsStr(
    Gene          = Gene,
    regionsFile   = regionsFile,
    regionsFormat = regionsFormat,
    geneCol       = geneCol,
    regionsCol    = regionsCol
  )

  # Obtain SCISSOR-style Ranges object for this gene
  Ranges <- SCISSOR::get_Ranges(
    Gene       = Gene,
    regions    = regions,
    outputType = "part_intron"
  )

  # Sanity checks for bamFilesList / caseIDList
  if (!all(studylist %in% names(bamFilesList))) {
    stop("All elements of `studylist` must be present in names(bamFilesList).")
  }
  if (!all(studylist %in% names(caseIDList))) {
    stop("All elements of `studylist` must be present in names(caseIDList).")
  }

  # Single-study mode
  if (length(studylist)==1L) {
    st <- studylist[1]

    res <- .construct_pileupStudy(
      Gene            = Gene,
      regions         = regions,
      Ranges          = Ranges,
      BAMfiles        = bamFilesList[[st]],
      caseIDs         = caseIDList[[st]],
      max_depth       = max_depth,
      strand.specific = strand.specific,
      nCores          = nCores
    )

    if (!is.null(outFile)) {
      outDir <- dirname(outFile)
      if (!dir.exists(outDir)) {
        dir.create(outDir, recursive=TRUE, showWarnings=FALSE)
      }

      pileup     <- res$pileup
      regions    <- res$regions
      Ranges     <- res$Ranges
      geneRanges <- Ranges
      save(pileup, regions, Ranges, geneRanges, file=outFile)
      return(invisible(outFile))
    }

    return(res)
  }

  # Multi-study mode
  res.list <- lapply(
    studylist,
    function(st) {
      .construct_pileupStudy(
        Gene            = Gene,
        regions         = regions,
        Ranges          = Ranges,
        BAMfiles        = bamFilesList[[st]],
        caseIDs         = caseIDList[[st]],
        max_depth       = max_depth,
        strand.specific = strand.specific,
        nCores          = nCores
      )
    }
  )
  names(res.list) <- studylist

  pileupList <- lapply(res.list, `[[`, "pileup")
  geneRanges <- Ranges

  if (!is.null(outFile)) {
    outDir <- dirname(outFile)
    if (!dir.exists(outDir)) {
      dir.create(outDir, recursive=TRUE, showWarnings=FALSE)
    }

    save(pileupList, regions, geneRanges, file=outFile)
    return(invisible(outFile))
  }

  list(pileupList=pileupList, regions=regions, geneRanges=geneRanges)
}


#' Core helper to build exon-only pileup
#'
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param cases a vector of specific samples among all samples in pileup. If NULL, all samples are selected. Default is NULL.
#' @param study a character of study abbreviation in the pileupList. Default is NULL.
#' @return a numeric matrix of exon-only coverage (rows: exon positions, columns: samples).
#' @references Choi, H.Y., Jo, H., Zhao, X. et al. SCISSOR: a framework for identifying structural changes in RNA transcripts. Nat Commun 12, 286 (2021).
#' @export

.build_pileupExon <- function(pileupPath, cases=NULL, study=NULL) {

  if (!file.exists(pileupPath)) {
    warning("File does not exist: ", pileupPath)
    return(NULL)
  }

  load(file=pileupPath)

  if (exists("pileup")) {
    # Case 1: file contains `pileup` directly
  } else if (exists("pileupList")) {
    # Case 2: file contains a list of pileups
    if (is.null(study)) {
      stop("`study` must be provided when the file contains `pileupList`.")
    }
    if (!study %in% names(pileupList)) {
      stop("Study '", study, "' not found in `pileupList`.")
    }
    pileup <- pileupList[[study]]
  } else {
    stop("Neither `pileup` nor `pileupList` found in loaded file: ", pileupPath)
  }

  if (!exists("regions")) {
    stop("`regions` object not found in loaded file; required for SCISSOR::build_pileup.")
  }

  # Keep exon location of union transcripts in pileup
  pileupData <- SCISSOR::build_pileup(
    Pileup     = pileup,
    regions    = regions,
    inputType  = "part_intron",
    outputType = "only_exon"
  )
  colnames(pileupData) <- colnames(pileup) # to keep the original sample IDs

  # Keep selected samples
  if (!is.null(cases)) {
    idx <- match(cases, colnames(pileupData))
    idx <- idx[!is.na(idx)]
    if (length(idx)==0L) {
      warning("None of the requested `cases` were found in column names.")
      return(pileupData[, 0, drop=FALSE]) # to return empty with correct row dimension
    }
    pileupData <- pileupData[, idx, drop=FALSE]
  }

  return(pileupData)
}


#' Get a focused pileup of exon location (for single-study)
#'
#' @param g the gene order in genelist
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param cases a vector of specific samples among all samples in pileup. If NULL, all samples are selected. Default is NULL.
#' @return a focused pileup is a the number of exon locations x the number of samples matrix for the g-th gene.
#' @export

get_pileupExon <- function(g, pileupPath, cases=NULL) {
  pileupData <- .build_pileupExon(
    pileupPath = pileupPath[g],
    cases      = cases,
    study      = NULL
  )
  return(pileupData)
}


#' Get a focused pileup of exon location (for multi-study)
#'
#' @param Gene a character of gene name
#' @param pileupPath file paths of coverage pileupData including .RData file names
#' @param cases a vector of specific samples among all samples in pileup. If NULL, all samples are selected. Default is NULL.
#' @param Study a character of study abbreviation in the pileupList. Default is NULL.
#' @return a focused pileup is a the number of exon locations x the number of samples matrix for the g-th gene or Gene.
#' @noRd

extract_pileupExon <- function(Gene, pileupPath, cases=NULL, Study=NULL) {
  pileupData <- .build_pileupExon(
    pileupPath = pileupPath,
    cases      = cases,
    study      = Study
  )
  return(pileupData)
}


#' Filter low expressed genes
#'
#' @param genelist a vector of gene names
#' @param TPM a gene expression counts matrix transformed by TPM
#' @param thr threshold. Default is 5.
#' @param pct percent. Default is 40.
#' @return a vector of filtered gene names
#' @examples
#' data("TOY_mrna")
#' filter_lowExpGenes(
#'   genelist = TOY_mrna$genes,
#'   TPM      = TOY_mrna$TPM
#' )
#' @export

filter_lowExpGenes <- function(genelist, TPM, thr=5, pct=40) {

  TPM2 <- stats::na.omit(TPM[match(genelist, rownames(TPM)), ])
  rows_to_keep <- apply(TPM2, 1, function(row) {mean(row<thr) < pct/100})
  genelist2 <- rownames(TPM2[rows_to_keep, , drop=FALSE])

  return(genelist2)
}
