#' A function to extract probe information from the human methylation and biomart dataset.
#'
#' This function use the output data from MRGN trio analysis and post filtering to extract additional information on probes for each trio. The extracted details are then used to calculate the difference between CpG in methylation probe and the midpoint of a nearby CpG island, CpG in methylation probe and the start of a gene, the gene length, methylation mean, GC content, and relation to CpG island.
#'
#' @param df A data frame generated from [postFilter]
#' @param modeltype A string with the name of causal model
#' @param TCGA.meth Methylation data matrix; probes in rows, individuals in columns
#' @param gene.exp Gene Expression data; genes in rows, individuals in columns
#' @param cna Copy Number Alteration data; genes in rows, individuals in columns
#' @param trios A trios data matrix with 4 columns
#' @param humanmeth The human methylation data matrix with probe ID, gene name, mapinfo, and cpg island information
#' @param biomart The biomart data matrix with gene name, gene start and gene end location, and GC content
#' @param nstartMeth The column where numeric values begin in methylation data
#' @param nstartGene The column where numeric values begin in gene expression data
#' @param nstartCNA The column where numeric values begin in copy number alteration data
#' @param type.ind A vector of positive or negative ER individuals ID
#'
#' @return A data frame with 19 columns:
#' @return IlmnID
#' @return      DNA Methylation Probe ID
#' @return UCSC_RefGene_Name
#' @return      Target gene name(s) where the Probe ID is found
#' @return Genome_Build
#' @return      The human reference genome information
#' @return CHR
#' @return      Chromosome containing the CpG
#' @return MAPINFO
#' @return      Chromosomal coordinates of the CpG
#' @return UCSC_CpG_Islands_Name
#' @return      Chromosomal coordinates of the CpG Island from UCSC.
#' @return Relation_to_UCSC_CpG_Island
#' @return      The location of the CpG relative to the CpG island.
#' @return trio_row
#' @return      The row number of the gene in the trio data
#' @return Biomart_GeneName
#' @return      Since a probe can be found in multiple genes, we look at one gene information at a time.
#' @return Gene.stable.ID
#' @return      The identifiers that refer to the same genomic features or species
#' @return Gene.stable.ID.version
#' @return      The increments when the set of transcripts linked to a gene changes
#' @return Gene_start
#' @return      The location where the gene begins
#' @return Gene_end
#' @return      The location where the gene ends
#' @return diff_cpG_mapinfo
#' @return      The difference between midpoint of nearby CpG location and the location of methylation probe
#' @return diff_mapinfo_geneStart
#' @return      The difference between location of methylation probe and the start of the gene
#' @return gene_length
#' @return      The length of the gene (difference between gene start and gene end)
#' @return Methyl_mean
#' @return      Mean of the methylation values for the corresponding probe ID
#' @return Methyl_sd
#' @return      Standard deviation of the methylation values for the corresponding probe ID
#' @return GC_content
#' @return      The content of Guanine-Cytosin for the probe ID.
#'
#' @export
#'
#' @seealso [findTrioAll] [addDupsCNA], [addDupsGENE] [findPCs] [findPCsGeneral] [analyzeTrios]
#'
#' @examples #Run postFilter() on results from analyzeTrios() to get an object 'df'
#' 
#' \dontrun{
#' probe_info_data = HumanMethProbeInfo(df, "M1.1", meth, gene, cna, 
#'         final.trios.df, humanmeth, biomart, 
#'         nstartMeth=5, nstartGene=3, nstartCNA=3, type.ind=clinical.pos)
#' }
#'


HumanMethProbeInfo <- function(df, modeltype, TCGA.meth, gene.exp, cna, trios, humanmeth, biomart, nstartMeth, nstartGene, nstartCNA, type.ind){
  
  #get the row number in trios data for model types
  trio_row <- df %>% filter(Inferred.Model3 == modeltype)%>% pull(trio_row)
  print(length(trio_row))
  #row numbers for those trios in the Methylation data
  meth_row_in_trios <- trios %>% slice(trio_row) %>% pull(meth.row)
  
  #finding common individuals among the 3 datasets
  com.ind <- colnames(gene.exp)[nstartGene:ncol(gene.exp)] %>%
    intersect(colnames(TCGA.meth)[nstartMeth:ncol(TCGA.meth)]) %>%
    intersect(colnames(cna)[nstartCNA:ncol(cna)])
  
  #finding common individuals between the 3 datasets and pos & neg ER individuals
  com.ind.type <- type.ind %>% unlist() %>% intersect(com.ind)
  
  #find the column number of the individuals
  ind.col.data.type = match(com.ind.type, colnames(TCGA.meth))
  
  #get the data from meth data using the row numbers
  
  dat <- TCGA.meth %>% slice(meth_row_in_trios)
  
  data1 <- dat %>% select(ind.col.data.type)
  
  #loop through each row in the data from meth data
  #initialize mean and standard deviation
  #loop through each row in the data from meth data
  #for(i in 1:nrow(data)){
  
  #get the mean for each row
  mean=rowMeans(data1,na.rm = T)
  #mean = mean(unlist(data[i,ind.col.data.type]), na.rm = TRUE)
  
  #save all the mean for all rows in data
  #final_mean = c(final_mean, mean)
  
  #get the standard dev for each row
  sd=apply(data1,1,sd,na.rm=T)
  #sd = sd(unlist(data[i,ind.col.data.type]), na.rm = TRUE)
  
  #save all the standard dev for all rows in data
  #final_sd = c(final_sd, sd)
  
  
  #final_mean <- apply(data[,ind.col.data.type],1,mean)
  #final_sd <- apply(data[,ind.col.data.type],1,sd)
  
  #match the probe id in data to human meth data
  rows_in_humanmeth<- match(dat$Row.names,humanmeth$Name)
  
  #extract the data from human meth using the row numbers
  final_humanmeth <- humanmeth %>% slice(rows_in_humanmeth)
  
  #final_data = cbind(final_humanmeth[,c("IlmnID","UCSC_RefGene_Name","Genome_Build", "CHR", "MAPINFO", "UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island")], trio_row)
  final_data = final_humanmeth %>% select("IlmnID","UCSC_RefGene_Name","Genome_Build", "CHR", "MAPINFO", "UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island") %>% mutate(Trio_row=trio_row)
  #remove the transcript and TSS columns for now
  uni_biomart <- biomart %>% distinct()
  #uni_biomart = unique(biomart[,-c(3,4,7)])
  #uni_biomart = biomart
  
  #initialize
  final = NULL
  final_res= NULL
  
  #loop through each row in the data
  for(i in 1:nrow(final_data)){
    print(i)
    avg_cpG_sites = NA
    diff_cpG_mapinfo = NA
    
    #match gene names between human meth and bio mart data
    genes_humanmeth <- unique(unlist(strsplit(as.character(final_data$UCSC_RefGene_Name[i]), ';')))
    
    #break the values in the cpG column
    cpG_sites = unlist(strsplit(as.character(final_data$UCSC_CpG_Islands_Name[i]), "[[:punct:]]"))
    
    if(length(cpG_sites) > 0){
      
      #find the average
      avg_cpG_sites = mean(c(as.numeric(cpG_sites[2]), as.numeric(cpG_sites[3])))
      
      ##### 1
      diff_cpG_mapinfo = avg_cpG_sites - final_data$MAPINFO[i]
      
    }
    
    for(j in 1:length(genes_humanmeth)){
      
      uni_biomart_row = which(uni_biomart$Gene.name == genes_humanmeth[j])
      
      if(length(uni_biomart_row) > 0){
        
        ##### 2
        diff_mapinfo_geneStart = final_data$MAPINFO[i] - uni_biomart$Gene.start..bp.[uni_biomart_row]
        
        #### 3
        gene_length = uni_biomart$Gene.end..bp.[uni_biomart_row] - uni_biomart$Gene.start..bp.[uni_biomart_row]
        
        res = cbind(genes_humanmeth[j], uni_biomart[uni_biomart_row,], diff_cpG_mapinfo, diff_mapinfo_geneStart, gene_length, mean[i], sd[i])
        colnames(res) = c("Biomart_GeneName", colnames(uni_biomart),"diff_cpG_mapinfo", "diff_mapinfo_geneStart", "gene_length","Methyl_mean", "Methyl_sd")
        
        final = rbind(final, cbind(final_data[i,], res))
      }
    }
    
  }
  final_res=final %>% filter(Biomart_GeneName %in% trios$Gene.name)
  
  return(final_res)
}
