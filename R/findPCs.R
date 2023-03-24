#' A function to calculate Principal Components (PCs) and significantly associated PCs with the data.
#'
#' This function uses Principal Component Analysis which can be used to analyze the high dimensional data while retaining the information that is important in the data. The PC score matrix is the coefficients of the linear combination of the initial variables that form the Principle Components and we find the PCs that are highly associated with the data.
#'
#' @param data A data matrix; probes or genes in rows, individuals in columns
#' @param startCol The column where numeric values begin in the data
#' @param GeneNameCol The column where the gene names are located
#' @param type.ind A vector of positive or negative ER individuals ID
#' @param com.ind A vector of common individuals between the datasets you want to find PCs for
#' @param type A string of individual type ("ER pos" or "ER neg")
#' @param bsize The number of columns of data to use in each block of correlation calculations
#'
#' @return A list of length 4:
#' @return pca.type
#' @return          A PC score matrix
#' @return data.with.conf
#' @return          A list of 7 elements with the first element being a list of significantly associated PCs to each gene or probe
#' @return new.data.no.var
#' @return          A updated data matrix used for analysis with individuals in rows and probes or genes in column
#'
#'
#' @export
#'
#' @seealso [prcomp()] used to calculate the PC score matrix; [get.conf.matrix()] used to calculate the significantly associated PCs
#'
#' @examples #load the datasets
#' @examples data(clinical.neg)
#' @examples data(clinical.pos)
#'
#' @examples #Find common individuals between the methylation and gene expression dataset
#' @examples com.ind = intersect(colnames(gene)[3:ncol(gene)], colnames(meth)[5:ncol(meth)])
#'
#' @examples #Use the function to get PC score matrix and significantly associated PCs.
#' @examples pc.gene = findPCs(gene, 3, 1, clinical.pos[,1], com.ind, "Pos", 1)
#'
#' @examples #The PC matrix
#' @examples pc.gene[[1]]
#'
#' @examples #List of significantly associated PCs
#' @examples pc.gene[[2]]$sig.asso.covs
#'
#' @import MRGN
#' @import usethis


findPCs <- function(data, startCol, GeneNameCol, type.ind, com.ind, type, bsize){

  #finding common individuals between the 3 datasets and pos & neg ER individuals
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
