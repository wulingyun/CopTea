#' Compute the node depths
#' 
#' Compute the node depths for given seed nodes and network
#' 
#' This function calculates the node depths, which is used in \code{\link{neeat}}.
#' 
#' @param seed.node.sets Logical matrix indicated the seed nodes associated with specific conditions (gene sets).
#' @param net The adjacent matrix of network.
#' @param max.depth Integer for the maximum depth that will be computed.
#' 
#' @return This function returns a matrix with the same dimensions as \code{seed.node.sets}, where the element
#' \code{[i, j]} represents the depth of node \code{i} under the condition (gene set) \code{j}.
#' 
#' @seealso \code{\link{neeat}}
#' 
#' @export
neeat_depths <- function(seed.node.sets, net, max.depth)
{
  if (is.null(dim(seed.node.sets)))
    seed.node.sets <- Matrix(as.logical(seed.node.sets))
  net.edges <- get_net_edges(net)
  fun <- function(i) get_node_depths(column(seed.node.sets, i), net.edges, max.depth)
  sapply(seq_len(dim(seed.node.sets)[2]), fun)
}


#' Compute the node depths with permutation
#' 
#' Compute the node depths for the given seed node set and its permutations.
#'
#' This function calculates the node depths of given seed nodes, togerther with the depths in random permutations of seed nodes.
#' 
#' @param seed.nodes Logical vector indicated the seed nodes.
#' @param net The adjacent matrix of network.
#' @param n.parm The number of permutations.
#' @param max.depth Integer for the maximum depth that will be computed.
#'
#' @return This function returns a matrix with the dimensions \code{c(length(seed.nodes), n.perm+1)}, where the element
#' \code{[i, 1]} represents the depth of node \code{i}, and \code{[i, j+1]} the depth in permutation \code{j}.
#' 
#' @seealso \code{\link{neeat}}, \code{\link{neeat_depths}}
#'
#' @export
neeat_depths_with_permutation <- function(seed.nodes, net, n.perm, max.depth = 1)
{
  n.all <- length(seed.nodes)
  n.seed <- sum(seed.nodes)
  vi <- integer(n.perm*n.seed)
  for (i in 1:n.perm)
    vi[(i-1)*n.seed + (1:n.seed)] <- sample.int(n.all, n.seed)
  vj <- rep(1:n.perm, each = n.seed)
  seed.perm <- sparseMatrix(vi, vj, x=T, dims=c(n.all, n.perm))
  neeat_depths(seed.perm, net, max.depth)
}


# Compute the node depths for a given seed node set
get_node_depths <- function(seed.nodes, net.edges, max.depth)
  .Call(NE_GetDepths, net.edges$edges, net.edges$index, seed.nodes, max.depth)


# Compute the edge depths from given node depths
get_edge_depth <- function(node.depth, net.edges, max.depth)
{
  edges <- net.edges$edges
  edges <- matrix(edges[edges[,1] < edges[,2], ], ncol=2)
  depth <- node.depth[edges[,1]] + node.depth[edges[,2]]
  depth[depth < 0] <- -1
  depth[depth > max.depth] <- -1
  depth
}


# Get the edge list
get_net_edges <- function(net, gene.filter = NULL)
{
  edges <- which(net != 0, arr.ind = T)
  if (!is.null(gene.filter))
    edges <- matrix(edges[gene.filter[edges[,1]] & gene.filter[edges[,2]], ], ncol=2)
  edges <- edges[order(edges[,1]), ]
  index <- findInterval(0:dim(net)[1], edges[,1])
  list(edges=edges, index=index)
}
