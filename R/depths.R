#' Compute the node depths
#' 
#' Compute the node depths for given core sets and network
#' 
#' This function calculates the node depths, which can be used for calling \code{\link{neeat}} in the
#' batch mode to reduce the computation time.
#' 
#' @param core.sets Logical matrix indicated the core genes associated with specific functions or pathways.
#' @param net The adjacent matrix of network.
#' @param max.depth Integer for the maximum depth that will be computed.
#' 
#' @return This function returns a matrix with the same dimensions as \code{core.sets}, where the element
#' \code{[i, j]} represents the depth of node \code{i} under the condition (function or pathway) \code{j}.
#' 
#' @seealso \code{\link{neeat}}
#' 
#' @export
neeat_depths <- function(core.sets, net, max.depth)
{
  net.edges <- net_edges(net)
  fun <- function(i) get_depths(column(core.sets, i), net.edges, max.depth)
  sapply(seq_len(dim(core.sets)[2]), fun)
}


get_depths <- function(core.set, net.edges, max.depth)
  .Call(NE_GetDepths, net.edges$edges, net.edges$index, core.set, max.depth)


edge_depth <- function(node.depth, net.edges, max.depth)
{
  edges <- net.edges$edges
  edges <- matrix(edges[edges[,1] < edges[,2], ], ncol=2)
  depth <- node.depth[edges[,1]] + node.depth[edges[,2]]
  depth[depth < 0] <- -1
  depth[depth > max.depth] <- -1
  depth
}
