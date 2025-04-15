#' LD Matrices for Top 10 Height-Associated Loci
#'
#' This dataset is a list of estimated linkage disequilibrium (LD) matrices, one for each of the top 10 loci associated with height.
#'
#' The loci were selected based on genome-wide significant associations (P < 5e-8) from a height GWAS. Clumping was performed using PLINK's clump-and-threshold (C+T) method with the following parameters:
#' \code{--clump-field P}, \code{--clump-kb 1000}, \code{--clump-p1 5e-8}, \code{--clump-p2 5e-8}, and \code{--clump-r2 0.9}, retaining independent index SNPs within each locus.
#'
#' LD matrices were estimated using a reference panel of 9,680 unrelated individuals (4,475 males and 5,205 females) randomly selected from the UK Biobank.
#'
#' Each element in the list corresponds to one locus and contains a symmetric numeric matrix of pairwise LD values (typically R-squared).
#'
#' @format A list of 10 numeric matrices. Each matrix is symmetric and corresponds to one genomic locus.
#' @usage data(LD)
#' @source UK Biobank; height GWAS summary statistics; PLINK clumping and LD estimation.
#' @keywords datasets LD height GWAS
#' @examples
#' data(LD)
#' length(LD)  # Number of loci
#' dim(LD[[1]])  # LD matrix for the first locus
"LD"
