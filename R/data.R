#' Toy mRNA-seq-like dataset for RNAshapeQC (SE input)
#'
#' A small synthetic \code{SummarizedExperiment} object mimicking mRNA-seq
#' coverage-based quality control (QC) inputs. It is used in the vignette to
#' demonstrate degradation-based metrics such as decay rate (DR), degradation
#' score (DS), and the degraded/intact index (DII).
#'
#' @format A \code{SummarizedExperiment} object with:
#' \describe{
#'   \item{assays}{Two matrices:
#'     \describe{
#'       \item{DR}{A numeric matrix of decay rates (genes x samples).}
#'       \item{TPM}{A numeric matrix of TPM expression values with the same
#'       dimension and dimnames as \code{DR}.}
#'     }}
#'   \item{rowData}{A \code{DataFrame} containing gene-level metadata,
#'     including \code{gene_length} (bp).}
#' }
#'
#' @details
#' The dataset contains 100 synthetic genes and 10 synthetic samples.
#' All values were generated solely for demonstration and testing purposes
#' and do not correspond to real biological data.
#'
#' @examples
#' data(TOY_mrna_se)
#' TOY_mrna_se
#'
#' @name TOY_mrna_se
#' @docType data
#' @keywords datasets
NULL


#' Toy total RNA-seq-like dataset for RNAshapeQC (SE input)
#'
#' A small synthetic \code{SummarizedExperiment} object mimicking total
#' RNA-seq coverage-based quality control (QC) inputs. It is used in the
#' vignette to demonstrate coverage-shape metrics such as mean coverage depth
#' (MCD), window coefficient of variation (wCV), and the suboptimal/optimal
#' index (SOI).
#'
#' @format A \code{SummarizedExperiment} object with:
#' \describe{
#'   \item{assays}{Two matrices:
#'     \describe{
#'       \item{MCD}{A numeric matrix of mean coverage depth (genes x samples).}
#'       \item{wCV}{A numeric matrix of window coefficients of variation with
#'       the same dimension and dimnames as \code{MCD}.}
#'     }}
#' }
#'
#' @details
#' The dataset contains 100 synthetic genes and 10 synthetic samples.
#' All values were generated solely for demonstration and testing purposes
#' and do not correspond to real biological data.
#'
#' @examples
#' data(TOY_total_se)
#' TOY_total_se
#'
#' @name TOY_total_se
#' @docType data
#' @keywords datasets
NULL


#' Toy mRNA-seq-like dataset for RNAshapeQC (matrix input)
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
#'   \item{genelength}{A gene length (bp) vector with names matching the row names of \code{DR}.}
#' }
#'
#' @details
#' All values are synthetic and were generated solely for demonstration and
#' testing. They do not correspond to any real samples or cohorts.
#'
#' @examples
#' data("TOY_mrna_mat")
#' str(TOY_mrna_mat)
#'
#' @name TOY_mrna_mat
#' @docType data
#' @keywords datasets
NULL


#' Toy total RNA-seq-like dataset for RNAshapeQC (matrix input)
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
#'
#' @details
#' All values are synthetic and were generated solely for demonstration and
#' testing. They do not correspond to any real samples or cohorts.
#'
#' @examples
#' data("TOY_total_mat")
#' str(TOY_total_mat)
#'
#' @name TOY_total_mat
#' @docType data
#' @keywords datasets
NULL
