#' A function to find rows that have duplicated Entrez IDs
#'
#' In genomic datasets, duplicated rows sometimes have either the same gene name and Entrez ID or the same values for all the individuals which introduces repetition in the data. Hence, we extract the duplicated rows so we can remove them later.
#'
#' @param data A data with genes in rows and individuals in columns
#'
#' @return a vector with row numbers
#' @export
#'
#' @examples #load the example datasets
#' @examples data(meth)
#' @examples data(cna)
#' @examples data(gene)
#'
#' @examples #use the function to get the row numbers in the dataset
#' @examples row = findDups(cna)
#' @examples row
#' @examples cna[row,1:5]
#'


# Finding duplicates in entrez id
findDups <- function(data){


  #gives rows for which there are no NA as entrez id
  nona.ent <- which(!is.na(data$Entrez_Gene_Id))

  #the entrez ids of the above rows
  nona.entid <- data$Entrez_Gene_Id[nona.ent]
  nona.entid[1:5]

  #checks if the entred ids have duplicates or not (TRUE/FALSE)
  d.ent <- duplicated(nona.entid)

  #find which rows have duplicates (gives only the first one not the appearances)
  dups.rows <- which(d.ent == TRUE)

  #pick the first duplicated row and find the entrez id
  entz.ids <- nona.entid[dups.rows]

  #initialize variable for duplicated rows in CNA data
  dups.row.data <- NULL

  #find the rows where this particular entrez id is duplicated
  for(i in 1:length(entz.ids)){

    #since there are multiple entrez ids, we loop through each one of them
    dups <- which(data$Entrez_Gene_Id == entz.ids[i])

    #save all the row numbers
    dups.row.data <- c(dups.row.data, dups)
  }

  #return the result
  return(dups.row.data)
}

