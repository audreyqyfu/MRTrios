#' A function to form trios between three variables of interest for multiple genes
#'
#' A function to find trios by matching methylation probes, gene expression, and copy number alteration (CNA) for a vector of genes.
#'
#' @param meth.data Methylation data; probes in rows, individuals in columns
#' @param cna.data Copy Number Alteration data; genes in rows, individuals in columns
#' @param gene.data Gene expression data; genes in rows, individuals in columns
#' @param nStartMeth The column where numeric values begin in methylation data
#' @param nStartGene The column where numeric values begin in gene expression data
#' @param nStartCNA The column where numeric values begin in CNA data
#'
#' @return A data matrix of four columns:
#' @return gene.name
#' @return           name of the unique gene
#' @return meth.row
#' @return           row number of the gene in methylation data
#' @return cna.row
#' @return           row number of the gene in cna data
#' @return gene.row
#' @return           row number of the gene in gene expression data
#' @export
#'
#' @seealso [findDups()] to find duplicated rows; [removeNA()] to find rows with all NA for numeric values; [trios()] to form trio for each gene using gene names, [entrezCNA()] to match trios in CNA data with Gene Expression data based on entrez IDs, [entrezGENE()] to match trios in Gene Expression data with CNA data based on entrez IDs
#'
#' @examples #use the function to form trios of genes between methylation, gene expression, and cna data
#' @examples trios.df = findTrioAll(meth, cna, gene, 5, 3, 3)
#' @examples trios.df
#'



findTrioAll <- function(meth.data, cna.data, gene.data, nStartMeth, nStartGene, nStartCNA) {

  #rows with duplicates in CNA data
  dup.CNA <- findDups(cna.data)

  #rows with duplicates in Gene Exp data
  dup.GENE <- findDups(gene.data)

  #rows with NAs for all individual in Methylation data
  na.METH <- removeNA(meth.data, nStartMeth)

  #rows with NAs for all individual in Gene Exp data
  na.GENE <- removeNA(gene.data, nStartGene)

  #rows with NAs for all individual in CNA data
  na.CNA <- removeNA(cna.data, nStartCNA)

  #since the data has multiple genes in one row separated by ";", we split them
  split.genes <- unlist(strsplit(as.character(meth.data$Gene_Symbol), ';'))

  #then find the unique genes from the methylation data
  uni <- na.omit(unique(split.genes))

  if(length(dup.CNA) > 0){

      #remove any empty values in the uni list and skip the duplicated rows in CNA data
      rows <- na.omit (match(cna.data$Hugo_Symbol[dup.CNA], uni))

      if(length(rows) > 0){

        uni <- uni[-na.omit(rows)]

      }

    }

  if(length(dup.GENE) > 0){

      #skip the duplicated rows in Gene Exp data
      rows <- na.omit (match(gene.data$Hugo_Symbol[dup.GENE], uni))

      if(length(rows) > 0){

        uni <- uni[-na.omit(rows)]

      }

    }

  unique.genes = uni

  #to find the row in the methylation data
  meth.genes <- strsplit(as.character(meth.data$Gene_Symbol), ';')

  #initialize the tmp variable
  tmp <- NULL

  #loop through each gene
  for(i in 1:length(unique.genes)){
    print(i)

    #apply the function with the provided data
    tmp <- rbind(tmp,trios(unique.genes[i], meth.data, cna.data, gene.data, meth.genes, na.METH, na.GENE, na.CNA))

  }

  # dim(tmp)
  #print(tmp[1:5,])

  #find the rows that have NAs for each data type
  meth.na.rows <- which(is.na(tmp[,2]))
  cna.na.rows <- which(is.na(tmp[,3]))
  gene.na.rows <- which(is.na(tmp[,4]))

  #loop through the rows that did not have a gene match in CNA data
  if (length (cna.na.rows) > 0) {
    for(i in 1:length(cna.na.rows)){
      
      #find the entrez id match and replace the NA
      tmp[cna.na.rows[i],3] <- entrezCNA(cna.na.rows[i], tmp, cna.data, gene.data)
    }
  }

  #print(tmp[1:5,])

  #loop through the rows that did not have a gene match in Gene Exp data
  if (length (gene.na.rows) > 0) {
    for(i in 1:length(gene.na.rows)){
      
      #find the entrez id match and replace the NA
      tmp[gene.na.rows[i],4] <- entrezGENE(gene.na.rows[i], tmp, cna.data, gene.data)
    }
  }

  #print(tmp[1:5,])

  #return the final data
  return(tmp)

}
