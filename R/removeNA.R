
#' A function to find rows that have all NA values in a data
#'
#' This function uses a data with genes in rows and individuals in columns to extract row numbers that have all NA values
#'
#' @param data A data with genes in rows and individuals in columns
#' @param nStart The column where numeric values begin
#'
#' @return a vector with row numbers
#' @export
#'
#' @examples #load the example datasets
#' @examples data(meth)
#'
#' @examples #use the function and get the row numbers in the dataset that have NAs
#' @examples removeNA(meth, 5)

#function to find rows that have all NA values in a data
removeNA <- function(data, nStart){

  #find the rows that have NA as values for all the individuals in Methylation data
  na.cols <- which(rowSums(is.na(data[,nStart:ncol(data)])) == ncol(data[,nStart:ncol(data)]))

  #return the result
  return(na.cols)

}
