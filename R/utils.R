#' Extract core sets information from GO annotation
#'
#' Generate the core sets matrix from an object of class "Go3AnnDbBimap" in Bioconductor annotation package.
#'
#' This function generates the \code{core.sets} matrix required by \code{\link{neeat}} method.
#' 
#' @param go.map An object of class "Go3AnnDbBimap".
#' @param evidence Vector of string to filter the GO annotations, could be "ALL" or a set of evidence codes.
#' @param category Vector of string to filter the GO categories, could be "ALL" or a set of GO categories.
#' @param gene.set Vector of string to filter the genes.
#' 
#' @return This function returns a sparse matrix as required by \code{\link{neeat}} method.
#' 
#' @seealso \code{\link{neeat}}
#'
#' @examples
#' 
#' \dontrun{
#' source("http://bioconductor.org/biocLite.R")
#' biocLite("org.Hs.eg.db")
#' library(org.Hs.eg.db)
#' x <- get_core_sets(org.Hs.egGO2ALLEGS)
#' }
#'
#' @import Matrix
#' 
#' @export
get_core_sets <- function(go.map, evidence = "ALL", category = "ALL", gene.set = NULL)
{
  term.table <- AnnotationDbi::toTable(go.map)
  term.table[,1] <- toupper(term.table[,1])
  
  if (evidence != "ALL")
    term.table <- term.table[term.table[,3] %in% evidence, ]
  
  if (category != "ALL")
    term.table <- term.table[term.table[,4] %in% category, ]
  
  if (!is.null(gene.set)) {
    all.gene <- gene.set    
  }
  else
    all.gene <- unique(term.table[,1])
  
  toMatrix(term.table, all.gene)
}


#' Transform the gene lists into the matrix of gene sets
#'
#' Generate the gene sets matrix from several gene lists.
#'
#' This function generates the \code{gene.sets} matrix required by \code{\link{neeat}} method.
#' The name of each gene set is given by the name of corresponding component of \code{gene.lists},
#' and assigned to the column of \code{gene.sets}.
#' 
#' @param gene.lists A list of vectors of gene ids.
#' @param all.gene A vector of all gene ids.
#' 
#' @return This function returns a sparse matrix as required by \code{\link{neeat}} method.
#' 
#' @seealso \code{\link{neeat}}
#'
#' @import Matrix
#' 
#' @export
get_gene_sets <- function(gene.lists, all.gene)
{
  gene.set <- Matrix(F, nrow=length(all.gene), ncol=length(gene.lists), dimnames=list(all.gene, names(gene.lists)))
  for (i in 1:length(gene.lists)) {
    genes <- as.character(gene.lists[[i]])
    genes <- genes[genes %in% all.gene]
    gene.set[genes, i] <- T
  }
  gene.set
}


#' Cohen's kappa score
#' 
#' Calculate Cohen's kappa score for two vectors.
#' 
#' This function calculate Cohen's kappa score for two logical vectors.
#' 
#' @param x1 The first logical vector
#' @param x2 The second logical vector
#' 
#' @return The Cohen's kappa score
#' 
#' @export
kappa_score <- function(x1, x2)
{
  if (length(x1) != length(x2)) stop("x1 and x2 are not of same length!")
  t <- length(x1)
  p1 <- sum(x1)
  p2 <- sum(x2)
  p <- sum(x1 == x2) / t
  e <- (p1*p2 + (t-p1)*(t-p2)) / t^2
  (p-e) / (1-e)
}
