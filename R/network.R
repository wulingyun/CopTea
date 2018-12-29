#' Build the gene similarity network from annoations
#' 
#' Build the gene similarity network from the core sets matrix by using Cohen's kappa statistic.
#' 
#' This function build the gene similarity network from the gene annotations, which is given
#' by a core sets matrix. The Cohen's kappa statistic is used to compute the similarity between
#' two genes.
#' 
#' @param core.sets Logical matrix indicated the core genes associated with specific functions or pathways.
#' The rows correspond to genes, while the columns represent the functions or pathways.
#' @param kappa.threshold The threshold to establish a link between two genes if their kappa score is larger
#' than the threshold.
#' 
#' @return This function returns a sparse adjacent matrix of the gene similarity network.
#' 
#' @import Matrix
#' 
#' @export
build_sim_net <- function(core.sets, kappa.threshold = 0.9)
{
  core.sets <- t(core.sets)
  n.gene <- dim(core.sets)[2]
  all.gene <- colnames(core.sets)
  net <- Matrix(F, nrow=n.gene, ncol=n.gene, dimnames=list(all.gene, all.gene))
  for (i in 1:(n.gene-1))
  {
    gene.i <- Corbi::column(core.sets, i)
    for (j in (i+1):n.gene)
    {
      gene.j <- Corbi::column(core.sets, j)
      k <- kappa_score(gene.i, gene.j)
      if (k >= kappa.threshold)
      {
        net[i, j] <- T
        net[j, i] <- T
      }
    }
  }
  net
}
