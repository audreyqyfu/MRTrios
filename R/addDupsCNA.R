#' A function to match trios with CNA data using the package "org.Hs.eg.db"
#'
#' After matching the trios based on gene names and entrez ID, we still have missing data in the trios data matrix so we match the genes with a package available with additional information on the gene.
#'
#' @param final A trios data matrix with 4 columns; gene name, methylation row, cna row, gene expression row
#' @param cna Copy Number Alteration dataset; genes in rows, individuals in columns
#'
#' @return a list with two data matrix
#' @return final
#' @return      The trios data matrix with NAs in the CNA column replaced with row number from the CNA dataset
#' @return dup.final
#' @return      A data matrix with trios for genes that are matched with multiple entrez IDs in the CNA data
#'
#' @export
#'
#' @seealso [findDups] to find duplicated rows
#'
#' @examples 
#' # Match trios using gene name and entrez ID
#' trios.df = findTrioAll(meth, cna, gene, 5, 3, 3)
#' trios.df
#'
#' @examples 
#' # Match additional entries in the CNA column of trios data using package "org.Hs.eg.db"
#' result = addDupsCNA(trios.df, cna)
#'
#' @examples 
#' # Initial trios data with entries in the CNA column filled in after matching
#' result[[1]]
#'
#' @examples 
#' # Additional trios data for genes that were matched 
#' # and had multiple entrez id matches in the CNA data
#' result[[2]]
#'


addDupsCNA <- function(final, cna){

  #initialize
  new.row <- NULL

  #rows with duplicates in CNA data
  dup.CNA <- findDups(cna)

  #which rows have NA in the cna column in trio file
  na.row.cna <- which(is.na(final[,3]) == TRUE)

  #load the r package
  xx <- as.list(org.Hs.egALIAS2EG)

  #find the unique genes that have NA in cna column
  genes.na.row <- unique(unlist(final[na.row.cna,1]))
  dup.final <- NULL


  for(i in 1:length(genes.na.row)){
    print(i)

    #find the gene name for the missing entrez id
    gene.match <- which(names(xx) == as.character(genes.na.row[i]))
    gene.match

    if(length(gene.match) > 0){

      #find the entrez id
      entz.id.new <- xx[[gene.match]]
      entz.id.new

      if(length(entz.id.new) == 1){

        #find the row number for that entrez id in cna data
        new.row <- which(cna$Entrez_Gene_Id == entz.id.new)
        new.row

      }else if(length(entz.id.new) > 1){
        new.row = NULL

        #since the length is greater than 1
        #we loop through each entrez id and find a row where the match is found
        for(j in 1:length(entz.id.new)){

          row <- which(cna$Entrez_Gene_Id == entz.id.new[j])
          new.row <- c(new.row, row)

        }
      }

      #if there is atleast one match found
      if(length(new.row) > 0){

        save.row <- NULL

        #we loop through each row and make sure they are not in duplicated row
        #if they are we skip over them
        for(k in 1:length(new.row)){

          if((new.row[k] %in% dup.CNA) == FALSE){
            save.row = c(save.row,new.row[k])
          }
        }

        new.row <- save.row
      }

      #check if length greater than 0
      if(length(new.row) == 1){
        print(i)

        #find which rows have the gene in cna data
        final.cna.row <- which(final[,1] == as.character(genes.na.row[i]))
        final[final.cna.row,3] <- new.row
      }

      #check if length greater than 1
      if(length(new.row) > 1){
        print(i)


        #find the rows in the cna data that contain the specified gene
        final.cna.row <- which(final[,1] == as.character(genes.na.row[i]))

        # we save all the remaining values (except the first one) to dup.final
        for(k in 2:length(new.row)){

          dup.final <- rbind(dup.final,final[final.cna.row,])
          dup.final[which(is.na(dup.final[,3]) == TRUE),3] <- new.row[k]

        }

        # and save the first value to the original file
        final[final.cna.row,3] <- new.row[1]

      }

    }
  }

  #return the two lists
  return(list(final,dup.final))

}
