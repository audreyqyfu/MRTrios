#' A general function to derive Principal Components (PCs) 
#' and identify significantly associated PCs with the input variables.
#'
#' This function performs principal component (PC) analysis on the input data
#' and identifies scores that are significantly associated with the variables in the input. 
#' 
#'
#' @param data A data matrix; probes or genes in rows, individuals in columns
#' @param startCol The column where numeric values begin in the data
#' @param GeneNameCol The column where the gene names are located
#' @param com.ind A vector of common individuals between the datasets you want to find PCs for
#' @param bsize The number of columns of data to use in each block of correlation calculations
#'
#' @return A list of length 3:
#' @return pca.type
#' @return          A PC score matrix.
#' @return data.with.conf
#' @return          A list of 6 elements:
#'                  sig.asso.covs (list of associated PCs for each gene), 
#'                  pvalues (data matrix), 
#'                  qvalues (list of one data frame), 
#'                  cors (data matrix), 
#'                  sig (list of one data frame containing TRUE/FALSE),
#'                  adj.p (data matrix).
#' @return new.data.no.var
#' @return          A updated data matrix with individuals in rows and probes/genes in columns.
#'
#'
#' @export
#'
#' @seealso [prcomp] used to calculate the PC score matrix; 
#'          [get.conf.matrix] used to identify PCs that are significantly associated with the input;
#'          [findPCs] used to derive and identify PCs for the breast cancer cohort with two subtypes (ER+ and ER-).
#'
#' @examples #load the datasets
#' @examples data(gene)
#' @examples data(meth)
#'
#' @examples #Find common individuals between the methylation and gene expression dataset
#' @examples com.ind = intersect(colnames(gene)[3:ncol(gene)], colnames(meth)[5:ncol(meth)])
#'
#' @examples #Use the function to get PC score matrix and significantly associated PCs.
#' @examples pc.gene = findPCsGeneral(as.data.frame(gene), 3, 1, com.ind, 1)
#'
#' @examples #The PC matrix
#' @examples dim (pc.gene[[1]])
#' @examples pc.gene[[1]]
#'
#' @examples #List of significantly associated PCs
#' @examples pc.gene[[2]]$sig.asso.covs
#'
#' @import MRGN
#' @importFrom stats complete.cases


#findPCs <- function(data, startCol, GeneNameCol, type.ind, com.ind, type, bsize){
findPCsGeneral <- function(data, startCol=3, GeneNameCol=1, com.ind, bsize){
    
  #finding common individuals between the 3 datasets and pos & neg ER individuals
#  com.ind.type <- intersect(unlist(type.ind), com.ind)

  #find the column number of the individuals
#  ind.col.data.type = match(com.ind.type, colnames(data))
  ind.col.data.type = match(com.ind, colnames(data))
  
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

  #find the rows in data that have atleast one NA values
  na.rows <- which(complete.cases(t(new.data.no.var)) == FALSE)

  ############################################################################################
  #
  #                   Get the PC scores
  #
  ##########################################################################################

  if(length(na.rows) > 0){

    #calculate the PC score matrix
    pca.type <- prcomp(new.data.no.var[,-na.rows], scale = TRUE)


  }else{

    pca.type <- prcomp(new.data.no.var, scale = TRUE)

  }

  dim(pca.type$x)

  #get the significant associated pcs
  data.with.conf = get.conf.matrix(new.data.no.var, pca.type$x, blocksize = bsize, apply.qval = TRUE, pi0.method = "bootstrap", measure = "correlation", adjust_by = "all")


  #replace the column numbers to corresponding gene names
  colnames(new.data.no.var) <- data.info
  dim(new.data.no.var)

  #return the 3 datasets
  return(list(pca.type$x, data.with.conf, new.data.no.var))

}
