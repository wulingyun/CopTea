#' Read network information from text file
#' 
#' Read the network information from a text file with specific format.
#' 
#' This function reads the network information from a text file with specific format:
#' each line contains two strings separated by spaces, which correspond to the 
#' names of two end points of one edge in the network.
#' 
#' @param file The name of text file
#' @return A list with the following components:
#'   \item{size}{The number of network nodes}
#'   \item{node}{The vector of network node names}
#'   \item{matrix}{The logical adjacency matrix}
#' 
#' @seealso \code{\link{write_net}}
#' 
#' @import Matrix
#' 
#' @export
read_net <- function(file)
{
  net.text <- as.matrix(read.table(file, fill=T, as.is=T, col.names=1:max(count.fields(file))))
  net.node <- unique(as.character(net.text))
  net.node <- net.node[net.node != ""]
  net.edge <- cbind(as.character(net.text[,1]), as.character(net.text[,-1]))
  net.edge <- net.edge[net.edge[,2] != "", ]
  net.size <- length(net.node)
  node.id <- seq_along(net.node)
  names(node.id) <- net.node
  net.matrix <- sparseMatrix(node.id[net.edge[,1]], node.id[net.edge[,2]], x=T, dims=c(net.size, net.size), dimnames=list(net.node, net.node))
  list(size=net.size, node=net.node, matrix=net.matrix)
}


#' Write network information to text file
#' 
#' Write the network information to a text file with specific format.
#' 
#' This function writes the network information to a text file with specific format:
#' each line contains two strings separated by spaces, which correspond to the
#' names of two end points of one edge in the network.
#' 
#' @param net A list as returned by \code{\link{read_net}}
#' @param file The name of text file
#' 
#' @seealso \code{\link{read_net}}
#' 
#' @import Matrix
#' 
#' @export
write_net <- function(net, file)
{
  net.edge <- which(net$matrix != 0, arr.ind=1)
  net.edge <- matrix(net$node[net.edge], ncol=2)
  write.table(net.edge, file, quote=F, row.names=F, col.names=F)
}


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
    gene.i <- column(core.sets, i)
    for (j in (i+1):n.gene)
    {
      gene.j <- column(core.sets, j)
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
