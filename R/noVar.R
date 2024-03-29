#' A function to find columns in a data that have no variance
#'
#' Data with no variance might be easy to predict but has no uncertainity so we do not include those columns in our analysis.
#'
#' @param data A dataset with numeric values; individuals in rows and probes or genes in column
#'
#' @return A vector of column numbers that either have 0 variance or 2 or less than 2 non-NA values
#' @export
#'
#' @examples 
#' # Transpose the data since individuals are in columns
#' data = t(meth[,5:ncol(meth)])
#'
#' # Returns column numbers
#' col = noVar(data)
#' col
#' data[,col]

noVar <- function(data){

  #find columns with no variance
  var.col <- apply(data, 2, var, na.rm = TRUE)

  #return index of cols with 0 variance
  na.var1 <- which(var.col == 0)

  t <- apply(is.na(data), 2, sum)

  na.var2 <- which(t >= (nrow(data)-2))

  na.var <- c(na.var1, na.var2)

  return(na.var)

}
