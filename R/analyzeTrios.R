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
#' @param nStartMeth The column where numeric values begin in methylation data
#' @param nStartGene The column where numeric values begin in gene expression data
#'
#' @return a dataframe of dimension 1 x 14 with the following columns:
#' @return b11
#' @return      the indicator 0,1 value for the conditional test T1 ~ V | T2,U

#' @return b12
#' @return      the indicator 0,1 value for the conditional test T1 ~ T2 | V,U

#' @return b21
#' @return      the indicator 0,1 value for the conditional test T2 ~ V | T1,U

#' @return b22
#' @return      the indicator 0,1 value for the conditional test T2 ~ T1 | V,U

#' @return V1:T1
#' @return      the indicator 0,1 value for the marginal test between V1 and T1

#' @return V1:T2
#' @return      the indicator 0,1 value for the marginal test between V1 and T2

#' @return pb11
#' @return      the p-value for the conditional test T1 ~ V | T2,U

#' @return pb12
#' @return      the p-value for the conditional test T1 ~ T2 | V,U

#' @return pb21
#' @return      the p-value for the conditional test T2 ~ V | T1,U

#' @return pb22
#' @return      the p-value for the conditional test T2 ~ T1 | V,U

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
#' @examples #generate the PC score matrix and significantly associated PCs in both methylation and gene expression data for ER+ individuals.
#' @examples #The number of methylation probes in the methylation dataset (included with the package) is small in this example. To generate realistic result, we need whole genome data for this function as well as the findPCs() function.
#' @examples \dontrun{ pc.meth = findPCs(meth, 5, 2, clinical.pos[,1], com.ind, "Pos", 1)
#'                     pc.gene = findPCs(gene, 3, 1, clinical.pos[,1], com.ind, "Pos", 1)
#'                     pc.meth.pos.tmp = pc.meth.pos[[1]]
#'                     pc.gene.pos.tmp = pc.gene.pos[[1]]
#'                     sig.pcs.meth.tmp = pc.meth.pos[[2]]$sig.asso.covs
#'                     sig.pcs.gene.tmp = pc.gene.pos[[2]]$sig.asso.covs
#'                     final.result = analyzeTrios(meth, gene, cna, final.trios.df, pc.meth.pos.tmp, pc.gene.pos.tmp, sig.pcs.meth.tmp, sig.pcs.gene.tmp, clinical.pos, meth.table.pos, gene.table.pos, 2, 3, 5, 3)
#'                     final.result}
#'
#' @examples #perform the MRGN trio inference with confounding variables
#' @examples #The PC data used here (pc.meth.pos, pc.gene.pos, sig.pcs.meth, and sig.pcs.gene) are generated from the whole genome data.
#' @examples final.result = analyzeTrios(meth, gene, cna, final.trios.df[1:5,], pc.meth.pos, pc.gene.pos, sig.pcs.meth, sig.pcs.gene, clinical.pos, meth.table.pos, gene.table.pos, age.col=2, race.col=3, nStartMeth=5, nStartGene=3)
#' @examples final.result


analyzeTrios <- function(TCGA.meth, gene.exp, cna, trios, pc.meth, pc.gene, meth.sig.asso.pcs, gene.sig.asso.pcs, clinical, meth.table, gene.table, age.col, race.col, nStartMeth, nStartGene){

  # find the common individuals between the 3 datasets
  # pc matrix has common individuals from meth and gene exp
  com.ind <- intersect(rownames(pc.meth),colnames(cna))

  #find the rows in clinical data for the common individuals
  rows.clinical <- match(com.ind, unlist(clinical[,1]))

  #extract the age and race for those individuals
  age <- clinical[rows.clinical,age.col]
  race <- clinical[rows.clinical,race.col]

  #find the rows for the common individuals in the resp datasets
  ind.col.cna <- match(com.ind, colnames(cna))
  ind.col.gene <- match(com.ind, colnames(gene.exp))
  ind.col.meth <- match(com.ind, colnames(TCGA.meth))

  #initialize tmp
  tmp <- NULL

  #begin the loop for rows in trios
  for(i in 1:nrow(trios)){

    #check if the values for cna and gene exp are NA or not
    if(is.na(trios[i,3]) == FALSE & is.na(trios[i,4]) == FALSE) {

      if(rowSums(is.na(gene.exp[as.numeric(trios[i,4]),nStartGene:ncol(gene.exp)])) != ncol(gene.exp[,nStartGene:ncol(gene.exp)])){

        #print the row number
        print(i)

        #extract data for each dataset and create the trio
        trio.cna = t(cna[as.numeric(trios[i,3]),ind.col.cna])
        trio.gene = t(gene.exp[as.numeric(trios[i,4]),ind.col.gene])
        trio.meth = t(TCGA.meth[as.numeric(trios[i,2]),ind.col.meth])

        #combine the matrix
        trio.mat = as.data.frame(cbind(as.numeric(trio.cna), as.numeric(trio.gene), as.numeric(trio.meth)))

        #extract the updated row number from the indices table
        row.sig.pcs.meth <- unlist(meth.table[as.numeric(trios[i,2]),2])
        row.sig.pcs.gene <- unlist(gene.table[as.numeric(trios[i,4]),2])

        #find the row numbers for the common individuals in the pc score matrix
        com.ind.meth <- match(com.ind, rownames(pc.meth))
        com.ind.gene <- match(com.ind, rownames(pc.meth))

        #extract the column numbers for sig pcs
        sig.pcs.cols.meth <- as.integer(unlist(meth.sig.asso.pcs[row.sig.pcs.meth]))
        sig.pcs.cols.gene <- as.integer(unlist(gene.sig.asso.pcs[row.sig.pcs.gene]))

        #remvove PCs greater than 50
        insig.meth <- which(sig.pcs.cols.meth > 50)

        if(length(insig.meth) > 0){

          sig.pcs.cols.meth <- sig.pcs.cols.meth[-insig.meth]

        }else{

          sig.pcs.cols.meth <- sig.pcs.cols.meth

        }


        insig.gene <- which(sig.pcs.cols.gene > 50)

        if(length(insig.gene) > 0){

          sig.pcs.cols.gene <- sig.pcs.cols.gene[-insig.gene]

        }else{

          sig.pcs.cols.gene <- sig.pcs.cols.gene

        }

        #extract the sig columns from the pc matrix with the common individuals
        sig.pc.gene <- pc.gene[com.ind.gene, sig.pcs.cols.gene]
        sig.pc.meth <- pc.meth[com.ind.meth, sig.pcs.cols.meth]

        #count the total number of sig pcs
        total.pc.count <- length(unlist(gene.sig.asso.pcs[row.sig.pcs.gene])) + length(unlist(meth.sig.asso.pcs[row.sig.pcs.meth]))

        #create matrix with the trios and the confounding variables
        final.mat <- cbind(trio.mat, sig.pc.gene, sig.pc.meth, age, race)

        #apply MRGN and infer the trio
        res = infer.trio(as.data.frame(final.mat), use.perm = TRUE, is.CNA = TRUE, nperms = 500)

        #combine the row number of trios, model type, and pc count
        final <- cbind(i, res, total.pc.count)

        #combine the value for each trio in rows
        tmp <- rbind(tmp, final)

      }
    }
  }

  #return the dataset
  return(tmp)

}
