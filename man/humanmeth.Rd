\name{humanmeth}
\alias{humanmeth}
\docType{data}
\title{
 The human methylation dataset
}
\description{
 A dataset containing information on each probe in the mehtylation dataset.
}
\usage{data("humanmeth")}
\format{
  A data frame with 6 observations on the following 33 variables.
  \describe{
    \item{\code{IlmnID}}{Unique identifier from the Illumina CG database.}
    \item{\code{Name}}{Unique identifier from the Illumina CG database.}
    \item{\code{AddressA_ID}}{The address ID for the probe specific for the A allele.}
    \item{\code{AlleleA_ProbeSeq}}{The sequence of the probe identified in the AddressA_ID column.}
    \item{\code{AddressB_ID}}{The address ID for the probe specific for the B allele.}
    \item{\code{AlleleB_ProbeSeq}}{The sequence of the probe identified in the AddressB_ID column.}
    \item{\code{Infinium_Design_Type}}{Infinium I (2 probes/locus) or Infinium II (1 probe/locus).}
    \item{\code{Next_Base}}{For Infinium I probes, the nucleotide immediately following the CpG. Blank for Infinium II.}
    \item{\code{Color_Channel}}{For Infinium I probes, the color channel of the Next_Base signal.}
    \item{\code{Forward_Sequence}}{Plus (+) strand (HapMap) sequence (5'-3') flanking the CG.}
    \item{\code{Genome_Build}}{The human reference genome information}
    \item{\code{CHR}}{Chromosome containing the CpG}
    \item{\code{MAPINFO}}{Chromosomal coordinates of the CpG}
    \item{\code{SourceSeq}}{The original, genomic sequence used for probe design prior to bisulfite conversion.}
    \item{\code{Chromosome_36}}{Chromosome containing the CpG}
    \item{\code{Coordinate_36}}{Chromosomal coordinates of the CpG}
    \item{\code{Strand}}{The Forward (F) or Reverse (R) designation of the Design Strand.}
    \item{\code{Probe_SNPs}}{rsid(s) of SNP(s) located 10–50 bases from the target CpG.}
    \item{\code{Probe_SNPs_10}}{rsid(s) of SNP(s) located ≤ 10 bases from the target CpG.}
    \item{\code{Random_Loci}}{CpG loci chosen randomly by consortium members during the design process are marked "True".}
    \item{\code{Methyl27_Loci}}{CpG’s carried over from the HumanMethylation27 array (95\% carryover) are marked “True”.}
    \item{\code{UCSC_RefGene_Name}}{Target gene names, from the UCSC database. Multiple listings of the same gene name indicate splice variants.}
    \item{\code{UCSC_RefGene_Accession}}{The UCSC accession numbers of the target transcripts. Accession numbers are in the same order as the target gene transcripts.}
    \item{\code{UCSC_RefGene_Group}}{Gene region feature category describing the CpG position, from UCSC. Features are listed in the same order as the target gene transcripts}
    \item{\code{UCSC_CpG_Islands_Name}}{Chromosomal coordinates of the CpG Island from UCSC.}
    \item{\code{Relation_to_UCSC_CpG_Island}}{The location of the CpG relative to the CpG island.}
    \item{\code{Phantom}}{Classifications from the Functional Annotation of the Mammalian Genome (FANTOM) consortium as a low- or high-CpG density region and chromosomal coordinates from the FANTOM consortium of enhancer regions associated with FANTOM5 promoters.}
    \item{\code{DMR}}{Differentially methylated regions}
    \item{\code{Enhancer}}{Predicted enhancer elements as annotated in the original 450K design are marked True.}
    \item{\code{HMM_Island}}{Hidden Markov Model Islands. Chromosomal map coordinates of computationally predicted CpG islands.}
    \item{\code{Regulatory_Feature_Name}}{Chromosomal map coordinates of the regulatory feature }
    \item{\code{Regulatory_Feature_Group}}{Description of the regulatory feature referenced in Regulatory_Feature_Name, as provided by the Methylation Consortium.}
    \item{\code{DHS}}{DNase I Hypersensitivity Site}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13534
 https://knowledge.illumina.com/microarray/general/microarray-general-reference_material-list/000001567
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(humanmeth)
}

\keyword{datasets}

\seealso{
\code{\link[=HumanMethProbeInfo]{HumanMethProbeInfo()}}
}


