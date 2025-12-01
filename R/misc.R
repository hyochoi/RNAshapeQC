#' Convert character columns to numeric columns in a data frame
#'
#' @param df a data frame
#' @param colnums the column numbers to be converted
#' @references https://www.geeksforgeeks.org/how-to-convert-dataframe-column-from-character-to-numeric-in-r/#
#' @noRd

convert_Chr2Numcol = function(df, colnums) {
  vec <- c(colnums)
  df[ , vec] <- apply(df[ , vec,drop=FALSE], 2, function(x) as.numeric(as.character(x)))

  return(df)
}


#' Convert a continuous variable to eCDF values from ecdf function
#'
#' @param data a continuous variable for x-axis in the eCDF plot; a 1 x the number of genes or samples matrix; rownames are geneSymbol or NewSampleId.
#' @param name a column name for eCDF Values (name_ecdf)
#' @param margin 1 and 2 return eCDF values for gene information and sample information, respectively.
#' @references https://statisticsglobe.com/extract-ecdf-values-from-function-r
#' @references https://stats.stackexchange.com/questions/30858/how-to-calculate-cumulative-distribution-in-r
#' @noRd

convert_Cont2eCDF = function(data, name, margin){

  # Create ecdf function
  fun_ecdf <- stats::ecdf(data)

  # Apply ecdf function
  my_ecdf <- fun_ecdf(data)

  # Combine x & eCDF values
  data_ecdf1 <- data.frame(data, my_ecdf)
  data_ecdf2 <- cbind(as.matrix(rownames(data_ecdf1)), data_ecdf1)
  data_ecdf3 <- data_ecdf2[,-c(2)]

  if (margin==1) {
    colnames(data_ecdf3) <- c("geneSymbol", sprintf("%s_ecdf",name))

  } else if (margin==2) {
    colnames(data_ecdf3) <- c("NewSampleId", sprintf("%s_ecdf",name))

  } else {
    stop(margin," is not an option for margin.")
  }

  return(data_ecdf3)
}


#' Convert a matrix to a data frame using pivot_longer
#'
#' @param mat a matrix with observations in rows, variables in columns, and values in cells
#' @param rcenames row, column, and cell names of a matrix
#' @references https://dcl-wrangle.stanford.edu/pivot-basic.html
#' @noRd

convert_pivot.longer = function(mat, rcenames) {

  row <- rownames(mat)
  mat1 <- data.frame(row, mat)
  mat2 <- tidyr::pivot_longer(mat1,
                              cols = !tidyr::starts_with("row"),
                              names_to = "col",
                              values_to = "cell")
  mat3 <- mat2[order(mat2$col, decreasing = FALSE), ]
  colnames(mat3) <- rcenames

  return(mat3)
}


#' Combine vectors as a matrix from objects
#'
#' @param filePath file paths including .RData file names
#' @param objName a character of object name
#' @param header logical; whether the text files have a header line.
#' @param skip integer; number of lines to skip before reading data from `.txt` files.
#' @param txtCol integer; column index in the text file that contains the
#'   numeric vector to be extracted.
#' @param margin 1 and 2 return for gene- and sample-level vectors, respectively.
#' @param rowNames a vector of gene names
#' @param colNames a vector of sample names
#' @param nCores the number of cores for parallel computing. Default is 32.
#' @return a gene x sample matrix
#' @examples
#' \dontrun{
#' DR.mat <- combine_vecObj(
#'   filePath = paste0("path/to/", genelist, "/", Study, "_", genelist, "_degradation.RData"),
#'   objName  = "allSampleDegRate",
#'   margin   = 1,
#'   rowNames = genelist,
#'   colNames = cases,
#'   nCores   = 32
#' )
#' }
#' @export

combine_vecObj = function(filePath, objName=NULL, header=NULL, skip=NULL, txtCol=NULL, margin, rowNames, colNames, nCores=32) {

  if (!(margin %in% c(1, 2))) {
    stop(margin, " is not an option for margin.")
  }

  # Make a list of vectors
  vecList <- parallel::mclapply(filePath, function(i) {
    ext <- tools::file_ext(i)

    if (ext=="RData") {
      if (is.null(objName)) stop("objName should be specified for RData files.")
      extract_RData(i, objName)

    } else if (ext=="txt") {
      if (is.null(header)||is.null(skip)||is.null(txtCol)) stop("header, skip, and txtCol should be specified for txt files.")
      df <- read.table(i, header=header, skip=skip)
      if (txtCol>ncol(df)) stop("txtCol exceeds number of columns in txt files: ", i)
      as.numeric(df[[txtCol]])

    } else {
      stop("Unsupported file type: ", ext)
    }
  }, mc.cores=nCores-1)

  # Combine by rows or columns
  if (margin==1) {
    mat <- do.call(rbind, vecList)
  } else if (margin==2) {
    mat <- do.call(cbind, vecList)
  }
  rownames(mat) <- rowNames
  colnames(mat) <- colNames

  return(mat)
}


#' Extract an object from .RData
#'
#' @param file .RData file
#' @param object object name
#' @return the object
#' @references https://stackoverflow.com/questions/65964064/programmatically-extract-an-object-from-collection-of-rdata-files
#' @noRd

extract_RData = function(file, object) {
  # Function for extracting an object from a .RData file created by R's save() command
  # Inputs: RData file, object name
  E <- new.env()
  load(file=file, envir=E)

  return(get(object, envir=E, inherits=F))
}


#' Repeat rows of a matrix
#'
#' @param mat a matrix
#' @param n the number of repeat rows
#' @references https://stackoverflow.com/questions/11121385/repeat-rows-of-a-data-frame
#' @noRd

repeach = function(mat, n){
  mat1 <- mat[rep(seq_len(nrow(mat)), each=n), ]
  mat2 = as.matrix(mat1)

  return(mat2)
}


#' Repeat copies of array (equivalent of repmat in Matlab)
#'
#' @param X a matrix
#' @param m m copies of X are arranged in row
#' @param n n copies of X are arranged in column
#' @references Hua Liu, Jinhong You & Jiguo Cao (2023). A Dynamic Interaction Semiparametric Function-on-Scalar Model, Journal of the American Statistical Association, 118:541, 360-373, DOI: 10.1080/01621459.2021.1933496
#' @noRd

repmat = function(X, m, n){
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X, mx, nx*n)), mx*m, nx*n, byrow=TRUE)
}


#' Adjust y-axis for a plot
#'
#' @references https://github.com/hyochoi/SCISSOR/blob/master/R/yaxis.hy.R
#' @noRd

yaxis.hy = function(mat){
  #  mat : d by n matrix
  tempmax <- max(mat) ;
  tempmin <- min(mat) ;
  templen <- tempmax-tempmin ;
  return(c(tempmin-0.002*templen, tempmax+0.002*templen)) ;
}
