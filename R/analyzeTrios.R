#' A function to infer the causal inference network with confounding variables using MRGN
#'
#' This function uses the trio data matrix with the confounding variables (Principal components and age and race of individuals). It returns the results of the coefficient and marginal test with the inferred model.
#'
#' @param TCGA.meth Methylation data matrix; probes in rows, individuals in columns
#' @param gene.exp Gene Expression data matrix; genes in rows, individuals in columns
#' @param cna Copy Number Alteration data matrix; genes in rows, individuals in columns
#' @param trios A trios data matrix with 4 columns; gene name, methylation row, cna row, gene expression row
#' @param pc.meth PC score matrix for the methylation data
#' @param pc.gene PC score matrix for the gene expression data
#' @param meth.sig.asso.pcs A list of significantly associated PCs with the methylation data for each trio
#' @param gene.sig.asso.pcs A list of significantly associated PCs with the gene expression data for each trio
#' @param clinical Clinical data with at least three columns: ID, age, and race
#' @param meth.table A indices table to track the changes made to the input methylation data while generating the PCs
#' @param gene.table A indices table to track the changes made to the input gene expression data while generating the PCs
#' @param age.col The column where the age is located in the clinical data
#' @param race.col The column where the race is located in the clinical data
#' @param nObs An integer specifying the minimum number of complete observations (no missing values) in a trio
#' @param nPCs An integer specifying the cutoff of the principal components to be selected.  For example, if nPCs=50, then the 50th PC and later ones in the score matrix are discarded.  Note that this is not the total number of PCs selected.
#' @param writeToFile TRUE if writing the trio analysis result to an external file while the function runs.  FALSE (default) otherwise.  If TRUE, need to specify the external file name.
#' @param file Name the file to write output to. 
#'
#' @return a dataframe of dimension 1 x 16 with the following columns:
#' @return index
#' @return      index of the trio
#' @return b11
#' @return      a binary value. 1 if the conditional test T1 ~ V | {T2,U} rejects the null, and 0 otherwise.

#' @return b12
#' @return      the binary outcome for the conditional test T1 ~ T2 | {V,U}

#' @return b21
#' @return      the binary outcome for the conditional test T2 ~ V | {T1,U}

#' @return b22
#' @return      the binary outcome for the conditional test T2 ~ T1 | {V,U}

#' @return V1:T1
#' @return      the binary outcome for the marginal test between V1 and T1

#' @return V1:T2
#' @return      the binary outcome for the marginal test between V1 and T2

#' @return pb11
#' @return      the p-value for the conditional test T1 ~ V | {T2,U}

#' @return pb12
#' @return      the p-value for the conditional test T1 ~ T2 | {V,U}

#' @return pb21
#' @return      the p-value for the conditional test T2 ~ V | {T1,U}

#' @return pb22
#' @return      the p-value for the conditional test T2 ~ T1 | {V,U}

#' @return pV1:T1
#' @return      the p-value for the marginal test between V1 and T1

#' @return pV1:T2
#' @return      the p-value for the marginal test between V1 and T2

#' @return Minor.freq
#' @return      the calculated frequency of the minor allele of the genetic variant
#'
#' @return Inferred.Model
#' @return      a string indicating the inferred model type as returned by class.vec()
#'
#' @return Total.PC.Count
#' @return      the number of PCs included in the trio analysis.
#' @export
#'
#' @seealso [infer.trio()] used to infer the causal models
#'
#' @examples #Find common individuals between the methylation and gene expression datasets
#' @examples data (cna)
#' @examples data (meth)
#' @examples data (gene)
#' @examples com.ind = intersect(colnames(gene)[3:ncol(gene)], colnames(meth)[5:ncol(meth)])
#'
#' @examples #match trios using gene name and entrez ID
#' @examples trios.df = findTrioAll(meth, cna, gene, 5, 3, 3)
#'
#' @examples #match additional entries in the CNA column of trios data using the package "org.Hs.eg.db"
#' @examples result = addDupsCNA(trios.df, cna)
#'
#' @examples #initial trios data with entries in the CNA column filled in after matching
#' @examples result[[1]]
#'
#' @examples #additional trios data for genes that were matched and had multiple entrez id matches in the CNA data
#' @examples result[[2]]
#'
#' @examples #use the function to match additional entries in the gene expression column of trios data using the package "org.Hs.eg.db"
#' @examples #It also merges the intial and additional trios data (res[[1]] and res[[2]]) and returns one data matrix
#' @examples final.trios.df = addDupsGENE(result[[1]], result[[2]], gene)
#' @examples final.trios.df
#'
#' @examples #generate the indices tables
#' @examples gene.table.pos = findIndex(gene, 3, 1, clinical.pos[,1], com.ind)
#' @examples meth.table.pos = findIndex(meth, 3, 1, clinical.pos[,1], com.ind)
#'
#' @examples #perform the MRGN trio inference with confounding variables
#' @examples #The PC data used here (pc.meth.pos, pc.gene.pos, sig.pcs.meth, and sig.pcs.gene) are pre-generated from the whole genome data.
#' @examples #See code at the end for an example of generating PC score matrices and associated PCs.
#' 
#' @examples final.result = analyzeTrios(meth, gene, cna, final.trios.df[1:5,], pc.meth.pos, pc.gene.pos, sig.pcs.meth, sig.pcs.gene, clinical.pos, meth.table.pos, gene.table.pos, age.col=2, race.col=3, nObs = 30, nPCs = 50)
#' @examples final.result
#' 
#' @examples #You can also write the results to an external file while the function runs
#' @examples 
#' 
#' \dontrun{ final.result = analyzeTrios(meth, gene, cna, final.trios.df[1:5,], pc.meth.pos, pc.gene.pos, sig.pcs.meth, sig.pcs.gene, clinical.pos, meth.table.pos, gene.table.pos, age.col=2, race.col=3, writeToFile = TRUE, file = "trio_results.txt")}
#'  
#' 
#' @examples #When the PC score matrix and significantly associated PCs are not available,
#' @examples #need to generate them first as follows. 
#' @examples #The number of methylation probes in the methylation dataset (included with the package) is small in this example. To generate realistic result, we need whole genome data for this function as well as the findPCs() function.
#' @examples 
#' 
#' \dontrun{ pc.meth = findPCs(meth, 5, 2, clinical.pos[,1], com.ind, "Pos", 1)
#'                     pc.gene = findPCs(gene, 3, 1, clinical.pos[,1], com.ind, "Pos", 1)
#'                     pc.meth.pos.tmp = pc.meth.pos[[1]]
#'                     pc.gene.pos.tmp = pc.gene.pos[[1]]
#'                     sig.pcs.meth.tmp = pc.meth.pos[[2]]$sig.asso.covs
#'                     sig.pcs.gene.tmp = pc.gene.pos[[2]]$sig.asso.covs
#'                     final.result = analyzeTrios(meth, gene, cna, final.trios.df, pc.meth.pos.tmp, pc.gene.pos.tmp, sig.pcs.meth.tmp, sig.pcs.gene.tmp, clinical.pos, meth.table.pos, gene.table.pos, age.col = 2, race.col = 3, nObs = 30, nPCs = 50)
#'                     final.result}
#'


analyzeTrios <- function(TCGA.meth, gene.exp, cna, trios, pc.meth, pc.gene, meth.sig.asso.pcs, gene.sig.asso.pcs, clinical, meth.table, gene.table, age.col, race.col, nObs = 30, nPCs = 50, writeToFile = FALSE, file){

  # find the common individuals between the 3 datasets
  # pc matrix has common individuals from meth and gene exp
  com.ind <- intersect(rownames(pc.meth),colnames(cna))

  #find the rows in clinical data for the common individuals
  rows.clinical <- match(com.ind, unlist(clinical[,1]))

  #extract the age and race for those individuals
  age <- as.data.frame(clinical)[rows.clinical,age.col]
  race <- as.data.frame(clinical)[rows.clinical,race.col]

  #find the rows for the common individuals in the resp datasets
  ind.col.cna <- match(com.ind, colnames(cna))
  ind.col.gene <- match(com.ind, colnames(gene.exp))
  ind.col.meth <- match(com.ind, colnames(TCGA.meth))
  
  # extract data for the common individuals
  cna.common <- as.data.frame(cna)[,ind.col.cna]
  gene.exp.common <- as.data.frame(gene.exp)[,ind.col.gene]
  meth.common <- as.data.frame(TCGA.meth)[,ind.col.meth]

  #find the row numbers for the common individuals in the pc score matrix
  com.ind.pc <- match(com.ind, rownames(pc.meth))

  #initialize result
  result <- NULL

  # write to file if writeToFile is TRUE
  if (writeToFile) {
    result.colnames <- c("Index", "b11", "b12", "b21", "b22", "V1:T1", "V1:T2", "pb11", "pb12", "pb21", "pb22", "pV1:T1", "pV1:T2", "Minor.freq", "Inferred.Model", "Total.PC.Count")
    write.table(t(result.colnames), file = file, sep = "\t", row.names = FALSE, col.names = FALSE, append = FALSE, quote = FALSE)
  }
  
  #begin the loop for rows in trios
  for(i in 1:nrow(trios)){

    #check if the values for cna and gene exp are NA or not
    if(!(is.na(trios[i,3]))  & !(is.na(trios[i,4]))) {

      #print the row number
      print(i)
      
      #extract data for each dataset and create the trio
      trio.cna = t(cna.common[as.numeric(trios[i,3]),])
      trio.gene = t(gene.exp.common[as.numeric(trios[i,4]),])
      trio.meth = t(meth.common[as.numeric(trios[i,2]),])
      
      #combine the matrix
      trio.mat = as.data.frame(cbind(as.numeric(trio.cna), as.numeric(trio.gene), as.numeric(trio.meth)))
      
      # check whether there are enough observations 
      na.counts <- rowSums(is.na (trio.mat))

      # since a row with NAs is ignored in regression,
      # need to make sure that there are enough observations left
      if(sum (na.counts == 0) > nObs){

        #extract the updated row number from the indices table
        row.sig.pcs.meth <- unlist(meth.table[as.numeric(trios[i,2]),2])
        row.sig.pcs.gene <- unlist(gene.table[as.numeric(trios[i,4]),2])

        #extract the column numbers for sig pcs
        sig.pcs.cols.meth <- as.integer(unlist(meth.sig.asso.pcs[row.sig.pcs.meth]))
        sig.pcs.cols.gene <- as.integer(unlist(gene.sig.asso.pcs[row.sig.pcs.gene]))

        #extract the sig columns from the pc matrix with the common individuals
        # also limit the PCs to those less than nPCs
        sig.pc.gene <- pc.gene[com.ind.pc, sig.pcs.cols.gene[which (sig.pcs.cols.gene < nPCs)]]
        sig.pc.meth <- pc.meth[com.ind.pc, sig.pcs.cols.meth[which (sig.pcs.cols.meth < nPCs)]]

        #count the total number of sig pcs
        Total.PC.Count <- length(unlist(gene.sig.asso.pcs[row.sig.pcs.gene])) + length(unlist(meth.sig.asso.pcs[row.sig.pcs.meth]))

        #create matrix with the trios and the confounding variables
        final.mat <- cbind(trio.mat, sig.pc.gene, sig.pc.meth, age, race)

        #apply MRGN and infer the trio
        res = infer.trio(as.data.frame(final.mat), use.perm = TRUE, is.CNA = TRUE, nperms = 500)

        #combine the row number of trios, model type, and pc count
        final <- cbind(i, res, Total.PC.Count)

        # write to file if writeToFile is TRUE
        if (writeToFile) {
          write.table(final, file = file, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
        }
        #combine the value for each trio in rows
        result <- rbind(result, final)

      }
    }
  }

  #return the dataset
  colnames(result)[1] <- "Index"
  return(result)

}
