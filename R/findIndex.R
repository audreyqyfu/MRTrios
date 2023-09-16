#' A function tracking index changes in data matrix after data cleaning
#'
#' After finding the individuals of interest in data, retain rows (probes or genes) of data that have less than three NAs.  Output the index of the row before and after this processing.
#'
#' @param data A data matrix with probes or genes in rows and individuals in columns.
#' @param startCol The column where numeric values begin in the data.
#' @param GeneNameCol The column where the gene names are located.
#' @param type.ind A vector of positive or negative ER individuals ID.  If there are no subtypes, use the same argument for 'com.ind' here.
#' @param com.ind A vector of common individuals shared across the datasets for which principal component analyses are performed.
#'
#' @return A data matrix with 2 columns; Column 1 is indices in the input data, Column 2 is the indices in data after removing columns with no variance
#'
#' @export
#'
#' @examples #load the datasets
#' @examples data(gene)
#' @examples data(meth)
#' @examples data(clinical.neg)
#' @examples data(clinical.pos)
#'
#' @examples #Find common individuals between the methylation and gene expression dataset
#' @examples com.ind = intersect(colnames(gene)[3:ncol(gene)], colnames(meth)[5:ncol(meth)])
#'
#' # Use the function to get the indices data matrix
#' gene.table.pos = findIndex(data = gene, startCol = 3, GeneNameCol = 1, 
#'         type.ind = clinical.pos[,1], com.ind = com.ind)
#'

findIndex <- function(data, startCol, GeneNameCol, type.ind, com.ind){

  #for breast cancer consisting of ER+ and ER- subtypes,
  #find common individuals between the individuals specified in type.ind
  #and the individuals shared across the genomic datasets
  com.ind.type <- intersect(unlist(type.ind), com.ind)

  #find the column number of the individuals
  ind.col.data.type = match(com.ind.type, colnames(data))

  #only save numeric values
  new.data.numeric <- t(data[,ind.col.data.type])

  gene.names.data <- data[,GeneNameCol]

  #get the columns with 0 variance and less than 3 non-NA values
  na.var.data <- noVar(new.data.numeric)

  if(length(na.var.data) > 0){

    #remove those columns from the data and the gene names
    data.info <- gene.names.data[-na.var.data]
    new.data.no.var <- new.data.numeric[,-na.var.data]

  }else{

    #if there are no columns with no variance, keep the data as it is
    data.info <- gene.names.data
    new.data.no.var <- new.data.numeric
  }



  ############################################################################################
  #
  #                   Create the indices tables
  #
  ##########################################################################################

  # The first column is the rows numbers in the original data
  # The second column is the updated row numbers in the updated data
  # after removing NA values and no variance columns
  col1 <- 1:nrow(data)
  col2 <- rep(NA, nrow(data))
  to.rem <- na.var.data

  if(length(to.rem) > 0){

    col2[-to.rem] <- 1:(nrow(data) - length(to.rem))

  }else{

    col2 = col1

  }

  col.mtx.data <- cbind(col1,col2)

  #return the 4 datasets
  return(col.mtx.data)

}
