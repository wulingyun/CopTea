#' Network enhanced enrichment analysis tool
#' 
#' Functional enrichment anslysis of gene sets by integrating network information
#' 
#' This is a novel functional enrichment analysis tool based on the gene set and network
#' information.
#' 
#' @param eval.gene.set Logical vector indicated the gene set for evaluating.
#' @param func.gene.sets Logical matrix indicated the functional gene sets associated with specific functions or pathways.
#' The rows correspond to genes, while the columns represent the functions or pathways.
#' @param net The adjacent matrix of network.
#' @param method A string indicated the NEEAT model, including "neeat", "neeat_hyper" and "hyper". 
#' Use "hyper" for traditional hypergeometric test, in which the network information is ignored.
#' @param max.depth Integer for the maximum depth considered in the NEEAT models.
#' @param rho The weight parameter for depths.
#' @param n.perm The number of permutations for calculating p-values. Minimum value 100 is required for the model "neeat".
#' Ignored in the models "neeat_hyper" and "hyper".
#' @param z.threshold The threshold for filtering out small Z-scores. The p-values are calculated only
#' for the Z-scores \code{>= z.threshold}, otherwise p-values are set as 1.
#' @param verbose Logical variable indicated whether output the full variables for computing p-value,
#' such as mean and variance.
#' @param n.cpu The number of CPUs/cores used in the parallel computation.
#' @param batch.size The desired size of batches in the parallel computation.
#' 
#' @return This function will return a 2-dimensional array of dimensions \code{c(5, dim(func.core.sets)[2])},
#' and each column \code{[,j]} containing the following components for the correponding gene set for evaluating
#' and functional gene set \code{func.gene.sets[,j]}:
#'   \item{\code{z.score}}{The Z-score for the gene set}
#'   \item{\code{p.value}}{The statistic significance for the gene set under specified NEEAT model}
#'   \item{\code{raw.score}}{The raw score for the gene set under specified NEEAT model}
#'   \item{\code{avg.score}}{The average score for random permutations of the gene set}
#'   \item{\code{var.score}}{The variance of scores for random permutations of the gene set}
#'
#' @seealso \code{\link{get_func_gene_sets}}, \code{\link{neeat_depths}}
#' 
#' @import Matrix parallel
#'
#' @export
neeat <- function(eval.gene.set, func.gene.sets, net = NULL,
                  method = "neeat", max.depth = 1, rho = 0.5, n.perm = 10000,
                  z.threshold = 0, verbose = FALSE, n.cpu = 1, batch.size = 5000)
{
  if (is.null(dim(func.gene.sets)))
    func.gene.sets <- Matrix(as.logical(func.gene.sets))
  n.perm <- max(100, n.perm)
  
  neeat.par <- new.env()
  neeat.par$max.depth = max.depth
  neeat.par$rho = rho
  neeat.par$n.perm = n.perm
  neeat.par$z.threshold = z.threshold
  neeat.par$verbose = verbose
  
  if (method == "neeat") {
    neeat.par$depths <- neeat_depths_with_permutation(eval.gene.set, net, n.perm, max.depth, n.cpu)
  }
  else if (method == "neeat_hyper")
    eval.gene.set[neeat_depths(eval.gene.set, net, max.depth, n.cpu) >= 0] <- T

  if (n.cpu > 1) {
    if (.Platform$OS.type == "windows")
      cl <- makeCluster(n.cpu)
    else
      cl <- makeForkCluster(n.cpu)
    n.cl <- length(cl)
    n.job <- dim(func.gene.sets)[2]
    n.batch <- max(1, round(n.job / (batch.size * n.cl))) * n.cl
    jobs <- splitIndices(n.job, n.batch)
    result <- clusterApply(cl, jobs, neeat_internal, eval.gene.set, func.gene.sets, net, method, neeat.par)
    stopCluster(cl)
    result <- unlist(result)
  }
  else {
    result <- neeat_internal(seq_len(dim(func.gene.sets)[2]), eval.gene.set, func.gene.sets, net, method, neeat.par)
  }
  if (neeat.par$verbose)
    output.names <- c("z.score", "p.value", "raw.score", "avg.score", "var.score")
  else
    output.names <- c("z.score", "p.value")
  result <- array(result, dim = c(length(output.names), dim(func.gene.sets)[2]),
                  dimnames <- list(output.names, colnames(func.gene.sets)))
  result
}

neeat_internal <- function(fgs.ids, gene.set, func.gene.sets, net, method, neeat.par)
{
  if (method == "neeat") {
    neeat_fgs <- function(i)
    {
      fgs <- column(func.gene.sets, i)
      neeat_score(fgs, neeat.par)
    }
  }
  else if (method == "hyper" || method == "neeat_hyper") {
    N <- dim(func.gene.sets)[1]
    M <- colSums(func.gene.sets)
    n <- sum(gene.set)
    neeat_fgs <- function(i)
    {
      m <- sum(column(func.gene.sets, i) & gene.set)
      neeat_hyper(N, n, M[i], m, neeat.par)
    }
  }
  else {
    stop("Incorrect parameters!")
  }
  sapply(fgs.ids, neeat_fgs)
}

neeat_score <- function(fgs, neeat.par)
{
  depths <- neeat.par$depths[fgs, , drop=F]
  score.matrix <- neeat.par$rho^depths
  score.matrix[depths < 0] <- 0
  scores <- colSums(score.matrix)
  raw.score <- scores[1]
  avg.score <- mean(scores[-1])
  var.score <- var(scores[-1])

  if (!is.nan(var.score) && var.score > 0)
    z.score <- (raw.score - avg.score) / sqrt(var.score)
  else
    z.score <- 0

  if (z.score >= neeat.par$z.threshold)
    p.value <- sum(scores[-1] >= raw.score) / neeat.par$n.perm
  else
    p.value <- 1.0

  if (neeat.par$verbose)
    c(z.score, p.value, raw.score, avg.score, var.score)
  else
    c(z.score, p.value)
}

neeat_hyper <- function(N, n, M, m, neeat.par)
{
  avg.score <- n * M / N
  var.score <- n * (M / N) * ((N - M) / N) * ((N - n) / (N - 1))
  
  if (var.score > 0)
    z.score <- (m - avg.score) / sqrt(var.score)
  else
    z.score <- 0
  
  p.value <- phyper(m-1, M, N-M, n, lower.tail=F)

  if (neeat.par$verbose)
    c(z.score, p.value, m, avg.score, var.score)
  else
    c(z.score, p.value)
}

