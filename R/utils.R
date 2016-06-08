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
