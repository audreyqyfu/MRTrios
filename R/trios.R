#' A function to form trios between three variables of interest for a gene
#'
#' A function to find trios by matching methylation probes, gene expression, and copy number alteration (CNA) for an individual gene.
#'
#' @param gene.name A unique gene name
#' @param meth Methylation data matrix; probes in rows, individuals in columns
#' @param cna Copy Number Alteration data; genes in rows, individuals in columns
#' @param gene Gene Expression data; genes in rows, individuals in columns
#' @param meth.genes List of all genes in the methylation data
#' @param na.meth a vector of rows in methylation data that are all NAs
#' @param na.gene a vector of rows in gene expression data that are all NAs
#' @param na.cna a vector of rows in CNA data that are all NAs
#'
#' @return a row vector of length 4:
#' @return gene.name
#' @return           name of the unique gene
#' @return meth.row
#' @return           row number of the gene in methylation data
#' @return cna.row
#' @return           row number of the gene in cna data
#' @return gene.row
#' @return           row number of the gene in gene expression data
#'
#' @export
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
#' @examples #Now we use the function with the required arguments
#'
#' @examples #Example 1
#'
#' @examples result = trios("C3orf35", meth, cna, gene, meth.genes, na.meth, na.gene, na.cna)
#' @examples result
#'
#' # Example 2
#'
#' # 'result' below should have NA in the gene.row column 
#' # because no match is found.
#' result = trios("OR4P4", meth, cna, gene, meth.genes, na.meth, na.gene, na.cna)
#' result
#'
#' # Example 3
#'
#' # 'result' below should to be NULL 
#' # since this gene has all NA values in the methylation dataset
#' result = trios("RBL2", meth, cna, gene, meth.genes, na.meth, na.gene, na.cna)
#' result
#'
#' @import splitstackshape
#' @import org.Hs.eg.db

#function to match trios (works for one gene at a time)
trios <- function (gene.name, meth, cna, gene, meth.genes, na.meth, na.gene, na.cna) {

  #assign row numbers for the list data
  #so the genes in the same list will have the same row number
  g1 <- rep(seq_along(meth.genes), sapply(meth.genes, length))
  #g1[1:5]

  #match the gene name and find the row in methylation data
  meth.row <- g1[which(unlist(meth.genes) == gene.name)]
  #meth.row

  #match the gene name and find the row in CNA data
  cna.row <- which(unlist(cna$Hugo_Symbol) == gene.name)
  #cna.row

  #if no match found assign NA
  if(length(cna.row) == 0){
    cna.row = NA
  }

  #match the gene name and find the row in Gene Exp data
  gene.row <- which(unlist(gene$Hugo_Symbol) == gene.name)
  gene.row

  #if no match found assign NA
  if(length(gene.row) == 0){
    gene.row = NA
  }

  #initialize results
  results = NULL

  #since there are multiple matches for each gene name in the meth data
  #we loop through each row and make sure it does not have NA for every individual
  for (i in 1:length(meth.row)){
    #print(i)
    for (j in 1:length (gene.row)) {
      for (k in 1:length (cna.row)) {
        #check for NA rows
        if((as.integer(meth.row[i]) %in% as.integer(na.meth)) == FALSE & (as.integer(gene.row[j]) %in% as.integer(na.gene)) == FALSE & (as.integer(cna.row[k]) %in% as.integer(na.cna)) == FALSE) {
          
          #print(i)
          #we save the gene name, i (which is the row number in meth data)
          #and row numbers in cna & gene data
          results <- rbind(results,as.data.frame(cbind(gene.name, as.numeric(meth.row[i]), as.numeric(cna.row[k]), as.numeric(gene.row[j]))))
        }
      }
    }


  }

  #if a result is found, assign column names
  if(length(results) > 0){

    colnames(results)[2] = "meth.row"
    colnames(results)[3] = "cna.row"
    colnames(results)[4] = "gene.row"

  }

  #return the result to the function
  return(results)
}
