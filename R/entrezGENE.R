#' A function to locate gene in Gene Expression data from CNA based on their entrez ID
#'
#' This function uses the trios matrix that we get from the trios() function and fills in the missing values by matching entrez ID in the Gene Expression data to the CNA data.
#'
#' @param gene.row row numbers in trio data matrix with NA for the Gene Expression row column
#' @param data a trios data matrix with 4 columns
#' @param cna Copy Number Alteration data; genes in rows, individuals in columns
#' @param gene Gene Expression data; genes in rows, individuals in columns
#'
#' @return a row number in the Gene Expression data to replace NA in the data
#'
#' @export
#'
#' @seealso [trios()] to get the input data (trios data) used in this function
#'
#' @examples #load the example datasets
#' @examples data(meth)
#' @examples data(cna)
#' @examples data(gene)
#'
#' @examples #Methylation data has multiple genes in one row separated by ";", so we split them
#' @examples meth.genes = strsplit(as.character(meth$Gene_Symbol), ';')
#'
#' @examples #get the rows with all NA values
#' @examples na.meth = removeNA(meth, 5)
#' @examples na.gene = removeNA(gene, 3)
#' @examples na.cna = removeNA(cna, 3)
#'
#' @examples #use the trios function to get the trio data
#' @examples result = trios("LRRC16A", meth, cna, gene, meth.genes, na.meth, na.gene, na.cna)
#' @examples result
#'
#' @examples #Extract which row has a missing value
#' @examples #Here we have one row in the data so it returns 1 but in case of larger data it will return a vector of row numbers
#' @examples row = as.integer(which(is.na(result[,4] == TRUE)))
#'
#' @examples #Now we use the function with the required arguments
#' @examples val = entrezCNA(row, result, cna, gene)
#' @examples val
#'
#' @examples #Then we replace the value in the trio data
#' @examples result[row,4] = entrezGENE(row, result, cna, gene)
#' @examples result
#'

#function to find entrez ID rows in Gene Exp data
entrezGENE <- function(gene.row, data, cna, gene){

  #use the gene row but take the column of the cna row
  entrez.row.gene <- which(gene$Entrez_Gene_Id == cna$Entrez_Gene_Id[as.integer(data[gene.row,3])])
  entrez.row.gene

  #check if the matched row is NA or has more than 1 match (duplicated)
  #If not proceed
  if(length(entrez.row.gene) != 1){

    #replace the matched rows with the respective NA spots in the tmp data
    entrez.row.gene <- NA
  }
  #return the data
  return(entrez.row.gene)

}
