#' A function to locate gene in CNA data from Gene Expression based on their entrez ID
#'
#' This function uses the trios matrix that we get from the trios() function and fills in the missing values by matching entrez ID in the CNA data to the Gene Expression data.
#'
#' @param cna.row row numbers in trio data matrix with NA for the CNA row column
#' @param data a trios data matrix with 4 columns
#' @param cna Copy Number Alteration data; genes in rows, individuals in columns
#' @param gene Gene Expression data; genes in rows, individuals in columns
#'
#' @return a row number in the CNA data to replace NA in the data
#'
#' @export
#'
#' @seealso [trios] to get the input data (trios data) used in this function
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
#' @examples result = trios("C13orf15", meth, cna, gene, meth.genes, na.meth, na.gene, na.cna)
#' @examples result
#'
#' @examples #Extract which row has a missing value
#' @examples #Here we have one row in the data so it returns 1 but in case of larger data it will return a vector of row numbers
#' @examples row = as.integer(which(is.na(result[,3] == TRUE)))
#'
#' @examples #Now we use the function with the required arguments
#' @examples val = entrezCNA(row, result, cna, gene)
#' @examples val
#'
#' @examples #Then we replace the value in the trio data
#' @examples result[row,3] = entrezCNA(row, result, cna, gene)
#' @examples result


#function to find entrez ID rows in CNA data
entrezCNA <- function(cna.row, data, cna, gene){

  #use the cna row but take the column of the gene row
  #column 4 has gene rows in tmp and column 3 has cna rows in tmp
  entrez.row.cna <- which(cna$Entrez_Gene_Id == gene$Entrez_Gene_Id[as.integer(data[cna.row,4])])
  entrez.row.cna

  #check if the matched row is NA (no match) or has more than 1 match (duplicated)
  #If not proceed
  if(length(entrez.row.cna) != 1){

    #replace the matched rows with the respective NA spots in the tmp data
    entrez.row.cna <- NA
  }
  #return the data
  return(entrez.row.cna)

}
