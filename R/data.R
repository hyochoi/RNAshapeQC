#' Toy mRNA-seq-like dataset for RNAshapeQC
#'
#' A small synthetic dataset mimicking mRNA-seq coverage-based quality control
#' (QC) inputs. It is used in the vignette to demonstrate degradation-based
#' metrics such as decay rate (DR), degradation score (DS), and the
#' degraded/intact index (DII).
#'
#' @format A list with 6 components:
#' \describe{
#'   \item{DR}{A numeric matrix of decay rates; genes in rows and samples in columns.
#'     In this toy dataset, it is a 100 (genes) x 10 (samples) matrix with
#'     row names like \code{"Gene001"} and column names like \code{"T01"}.}
#'   \item{genes}{A character vector of length 100 containing gene IDs used as
#'     row names in \code{DR} and \code{TPM}.}
#'   \item{samples}{A character vector of length 10 containing sample IDs used
#'     as column names in \code{DR} and \code{TPM}.}
#'   \item{protocol}{A single character string indicating the protocol used,
#'     here \code{"mRNA-seq"}.}
#'   \item{TPM}{A numeric matrix of TPM values; same dimension and dimnames
#'     as \code{DR}.}
#'   \item{genelength.mat}{A one-column numeric matrix of gene lengths (bp),
#'     with row names matching the row names of \code{DR}.}
#' }
#' @return
#' A named list containing synthetic mRNA-seq-like inputs for RNAshapeQC.
#'
#' @details
#' All values are synthetic and were generated solely for demonstration and
#' testing. They do not correspond to any real samples or cohorts.
#'
#' @usage data("TOY_mrna")
#'
#' @examples
#' data("TOY_mrna")
#' str(TOY_mrna)
#'
#' @docType data
#' @keywords datasets
"TOY_mrna"


#' Toy total RNA-seq-like dataset for RNAshapeQC
#'
#' A small synthetic dataset mimicking total RNA-seq coverage-based quality
#' control (QC) inputs. It is used in the vignette to demonstrate coverage-shape
#' metrics such as mean coverage depth (MCD), window coefficient of variation
#' (wCV), and the suboptimal/optimal index (SOI).
#'
#' @format A list with 5 components:
#' \describe{
#'   \item{MCD}{A numeric matrix of mean coverage depth; genes in rows and
#'     samples in columns. In this toy dataset, it is a 100 (genes) x 10
#'     (samples) matrix with row names like \code{"Gene001"} and column names
#'     like \code{"A01"}.}
#'   \item{wCV}{A numeric matrix of window coefficients of variation; same
#'     dimension and dimnames as \code{MCD}.}
#'   \item{genes}{A character vector of length 100 containing gene IDs used as
#'     row names in \code{MCD} and \code{wCV}.}
#'   \item{samples}{A character vector of length 10 containing sample IDs used
#'     as column names in \code{MCD} and \code{wCV}.}
#'   \item{protocol}{A single character string indicating the protocol used,
#'     here \code{"total RNA-seq"}.}
#' }
#' @return
#' A named list containing synthetic total RNA-seq-like inputs for RNAshapeQC.
#'
#' @details
#' All values are synthetic and were generated solely for demonstration and
#' testing. They do not correspond to any real samples or cohorts.
#'
#' @usage data("TOY_total")
#'
#' @examples
#' data("TOY_total")
#' str(TOY_total)
#'
#' @docType data
#' @keywords datasets
"TOY_total"
