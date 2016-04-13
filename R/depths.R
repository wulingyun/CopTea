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
#' @param keep.degree Logical variable indicated whether keeping the similar degree distribution in permutaions.
#'
#' @return This function returns a matrix with the dimensions \code{c(length(seed.nodes), n.perm+1)}, where the element
#' \code{[i, 1]} represents the depth of node \code{i}, and \code{[i, j+1]} the depth in permutation \code{j}.
#' 
#' @seealso \code{\link{neeat}}, \code{\link{neeat_depths}}
#'
#' @export
neeat_depths_with_permutation <- function(seed.nodes, net, n.perm, max.depth = 1, keep.degree = T)
{
  n.all <- length(seed.nodes)
  n.seed <- sum(seed.nodes)
  if (keep.degree) {
    k.nb <- min(10, n.all-n.seed)
    d.all <- colSums(net)
    d.seed <- d.all[seed.nodes]
    d.other <- d.all[!seed.nodes]
    perm.list <- lapply(d.seed, get_perm_candidates, d.all, d.other, k.nb)
    vi <- replicate(n.perm, sample_list_without_replace(perm.list))
  }
  else {
    vi <- replicate(n.perm, sample.int(n.all, n.seed))
  }
  vj <- rep(1:n.perm, each = n.seed)
  seed.perm <- sparseMatrix(vi, vj, x=T, dims=c(n.all, n.perm))
  neeat_depths(seed.perm, net, max.depth)
}

get_perm_candidates <- function(d, d.all, d.other, k.nb)
{
  d.diff <- abs(d.other - d)
  d.range <- d.diff[order(d.diff)][k.nb]
  d.max <- d + d.range
  d.min <- d - d.range
  which(d.all >= d.min & d.all <= d.max)
}

sample_list_without_replace <- function(perm.list)
{
  perm <- sapply(perm.list, sample_i, 1)
  dup <- duplicated(perm)
  while (sum(dup) > 0) {
    print(perm)
    p.dup <- perm.list[dup]
    p.dup <- lapply(p.dup, setdiff, perm)
    perm[dup] <- sapply(p.dup, sample_i, 1)
    dup <- duplicated(perm)
  }
  perm
}

sample_i <- function(x, ...)
{
  if (length(x) > 0)
    x[sample.int(length(x), ...)]
  else
    stop("Attempt to sample from empty set!")
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
