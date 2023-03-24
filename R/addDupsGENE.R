#' A function to match trios with Gene Expression data using the package "org.Hs.eg.db"
#'
#' After matching the trios based on gene names and entrez ID, we still have missing data in the trios data matrix so we match the genes with a package available with additional information on the gene.
#'
#' @param final A trios data matrix with 4 columns; gene name, methylation row, cna row, gene expression row
#' @param dup.final A data matrix with trios for genes that are matched with multiple entrez IDs; gene name, methylation row, cna row, gene expression row
#' @param cna Gene Expression dataset; genes in rows, individuals in columns
#'
#' @return a list with two data matrix
#' @return final.res
#' @return      The trios data matrix with both final and dup.final data merged with NAs replaced by row numbers in Gene Expression column
#'
#' @export
#'
#' @seealso [findDups()] to find duplicated rows
#'
#' @examples #use the function to match trios using gene name and entrez ID
#' @examples trios.df = findTrioAll(meth, cna, gene, 5, 3, 3)
#' @examples trios.df
#'
#' @examples #use the function to match additional entries in the CNA column of trios data using the package "org.Hs.eg.db"
#' @examples result = addDupsCNA(trios.df, cna)
#'
#' @examples #initial trios data with entries in the CNA column filled in after matching
#' @examples result[[1]]
#'
#' @examples #additional trios data for genes that were matched and had multiple entrez id matches in the CNA data
#' @examples result[[2]]
#'
#' @examples #use the function to match additional entries in the Gene Expression column of trios data using the package "org.Hs.eg.db"
#' @examples #It also merges the intial and additional trios data (res[[1]] and res[[2]]) and returns one data matrix
#' @examples final.trios.df = addDupsGENE(result[[1]], result[[2]], gene)
#' @examples final.trios.df
#'


addDupsGENE <- function(final, dup.final, gene){

  #initialize
  new.row <- NULL

  #rows with duplicates in Gene Exp data
  dup.GENE <- findDups(gene)

  #which rows have NA in the gene exp column in trio file
  na.row.gene <- which(is.na(final[,4]) == TRUE)

  #load the r package
  xx <- as.list(org.Hs.egALIAS2EG)

  #find the unique genes that have NA in gene exp column
  genes.na.row <- unique(unlist(final[na.row.gene,1]))

  for(i in 1:length(genes.na.row)){
    print(i)

    #find the gene name for the missing entrez id
    gene.match <- which(names(xx) == as.character(genes.na.row[i]))
    gene.match

    #check if there is match in the R library for the specified gene name
    if(length(gene.match) > 0){

      #find the entrez id
      entz.id.new <- xx[[gene.match]]
      entz.id.new

      #check the length of the entz id and find the corresponding row in the data
      if(length(entz.id.new) == 1){

        #find the row number for that entrez id in cna data
        new.row <- which(gene$Entrez_Gene_Id == entz.id.new)
        new.row

        #check if length is greater than 1
      }else if(length(entz.id.new) > 1){
        new.row = NULL

        #loop through the entz id to find the match in gene exp data
        for(j in 1:length(entz.id.new)){

          #find the match for each entz id
          row <- which(gene$Entrez_Gene_Id == entz.id.new[j])

          #save every row to new.row
          new.row <- c(new.row, row)

        }
      }

      #if there is atleast one match found
      if(length(new.row) > 0){

        save.row <- NULL

        #we loop through each row and make sure they are not in duplicated row
        #if they are we skip over them
        for(k in 1:length(new.row)){

          if((new.row[k] %in% dup.GENE) == FALSE){
            save.row = c(save.row,new.row[k])
          }
        }

        new.row <- save.row
      }

      #check if length greater than 0
      if(length(new.row) == 1){


        #find which rows have the gene in gene.exp data
        final.gene.row <- which(final[,1] == as.character(genes.na.row[i]))
        final[final.gene.row,4] <- new.row

      }

      #check if length greater than 1
      if(length(new.row) > 1){

        #find the rows in the gene exp data that contain the specified gene
        final.gene.row <- which(final[,1] == as.character(genes.na.row[i]))

        #finding the rows for the gene in dup.final since we already created this data
        # while finding rows for cna data
        dup.gene.row <- which(dup.final[,1] == as.character(genes.na.row[i]))

        #length decides if the gene is already in dup.final or not
        # if the length is greater than 0, we replace NA for gene.row in the existing data
        if(length(dup.gene.row) > 0){

          # save the first value to the original file
          final[final.gene.row,4] <- new.row[1]

          #for each new.row value starting from 2 save to dup.final
          for(k in 2:length(new.row)){

            dup.final[dup.gene.row[1:length(final.gene.row)],4] <- new.row[k]

          }
          # if it is not we create dup.final for the gene
        }else{

          #for each new.row value starting from 2 save to dup.final
          for(k in 2:length(new.row)){

            #copy the corresponding data from final to dup.final
            dup.final <- rbind(dup.final,final[final.gene.row,])

            #replace the NA with the new.row in gene.row column
            dup.final[which(is.na(dup.final[,4]) == TRUE),4] <- new.row[k]

          }

          # save the first value to the original file
          final[final.gene.row,4] <- new.row[1]

        }
      }
    }
  }

  #combine the datasets into one
  final.res = rbind(final,dup.final)

  #rename the second column
  colnames(final.res)[2] = "meth.row"

  #return the merged dataset
  return(final.res)


}
