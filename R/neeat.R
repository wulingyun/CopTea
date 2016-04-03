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
#' @param perm.batch The desired size of permutation batches in the parallel computation.
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
neeat <- function(eval.gene.set, func.gene.sets, net,
                  method = "neeat", max.depth = 1, rho = 0.5, n.perm = 100000,
                  verbose = FALSE, n.cpu = 1, perm.batch = 10000)
{
  if (!is.null(dim(eval.gene.set)))
    eval.gene.set <- eval.gene.set[, 1]
  if (is.null(dim(func.gene.sets)))
    func.gene.sets <- Matrix(as.logical(func.gene.sets))
  if (length(eval.gene.set) != dim(func.gene.sets)[1])
    stop("length(eval.gene.set) != dim(func.gene.sets)[1]!")
  max.depth <- max(0, max.depth)
  n.perm <- max(100, n.perm)
  
  neeat.par <- new.env()
  neeat.par$max.depth = max.depth
  neeat.par$rho = rho
  neeat.par$n.perm = n.perm
  neeat.par$verbose = verbose
  neeat.par$n.cpu = n.cpu
  neeat.par$perm.batch = perm.batch
  
  if (method == "neeat") {
    neeat.par$rho <- c(0, rho^(0:max.depth))
    depths <- neeat_depths(eval.gene.set, net, max.depth)
    raw.score <- as.numeric(neeat.par$rho[depths + 2] %*% func.gene.sets)
    if (n.cpu > 1) {
      if (.Platform$OS.type == "windows")
        cl <- makeCluster(n.cpu)
      else
        cl <- makeForkCluster(n.cpu)
      n.cpu <- length(cl)
      neeat.par$n.perm <- ceiling(n.perm / n.cpu)
      perm.score <- clusterCall(cl, neeat_perm, eval.gene.set, func.gene.sets, net, raw.score, neeat.par)
      stopCluster(cl)
      perm.score <- matrix(rowMeans(matrix(unlist(perm.score), ncol=n.cpu)), nrow=3)
    }
    else {
      perm.score <- neeat_perm(eval.gene.set, func.gene.sets, net, raw.score, neeat.par)
    }
    z.score <- (raw.score - perm.score[2, ]) / sqrt(perm.score[3, ])
    if (neeat.par$verbose)
      result <- rbind(z.score, perm.score, raw.score)
    else
      result <- rbind(z.score, perm.score[1, ])
  }
  else if (method == "neeat_hyper") {
    eval.gene.set[neeat_depths(eval.gene.set, net, max.depth) >= 0] <- T
    result <- neeat_hyper(eval.gene.set, func.gene.sets, neeat.par)
  }
  else if (method == "hyper") {
    result <- neeat_hyper(eval.gene.set, func.gene.sets, neeat.par)
  }
  else {
    stop("Incorrect method!")
  }
  
  if (neeat.par$verbose)
    output.names <- c("z.score", "p.value", "avg.score", "var.score", "raw.score")
  else
    output.names <- c("z.score", "p.value")
  result <- array(result, dim = c(length(output.names), dim(func.gene.sets)[2]),
                  dimnames <- list(output.names, colnames(func.gene.sets)))
  result
}


neeat_perm <- function(eval.gene.set, func.gene.sets, net, raw.score, neeat.par)
{
  n.batch <- ceiling(neeat.par$n.perm / neeat.par$perm.batch)
  n.perm <- ceiling(neeat.par$n.perm / n.batch)
  perm.score <- list()
  for (i in 1:n.batch) {
    depths <- neeat_depths_with_permutation(eval.gene.set, net, n.perm, neeat.par$max.depth)
    score.matrix <- matrix(neeat.par$rho[depths + 2], length(eval.gene.set), n.perm)
    perm.score[[i]] <- sapply(1:dim(func.gene.sets)[2], neeat_i, score.matrix, func.gene.sets, raw.score, n.perm)
  }
  matrix(rowMeans(matrix(unlist(perm.score), ncol=n.batch)), nrow=3)
}


neeat_i <- function(i, score.matrix, func.gene.sets, raw.score, n.perm)
{
  perm.score <- .Call(NE_ColSums, score.matrix, which(column(func.gene.sets, i)))
  avg.score <- mean(perm.score)
  var.score <- var(perm.score)
  p.value <- sum(perm.score >= raw.score[i]) / n.perm
  c(p.value, avg.score, var.score)
}


neeat_hyper <- function(gene.set, func.gene.sets, neeat.par)
{
  N <- dim(func.gene.sets)[1]
  M <- colSums(func.gene.sets)
  n <- sum(gene.set)
  jobs <- seq_len(dim(func.gene.sets)[2])
  if (neeat.par$n.cpu > 1) {
    if (.Platform$OS.type == "windows")
      cl <- makeCluster(neeat.par$n.cpu)
    else
      cl <- makeForkCluster(neeat.par$n.cpu)
    result <- parSapply(cl, jobs, neeat_hyper_i, N, n, M, gene.set, func.gene.sets, neeat.par)
    stopCluster(cl)
  }
  else {
    result <- sapply(jobs, neeat_hyper_i, N, n, M, gene.set, func.gene.sets, neeat.par)
  }
  result
}


neeat_hyper_i <- function(i, N, n, M, gene.set, func.gene.sets, neeat.par)
{
  M <- M[i]
  m <- sum(column(func.gene.sets, i) & gene.set)

  avg.score <- n * M / N
  var.score <- n * (M / N) * ((N - M) / N) * ((N - n) / (N - 1))
  
  if (var.score > 0)
    z.score <- (m - avg.score) / sqrt(var.score)
  else
    z.score <- 0
  
  p.value <- phyper(m-1, M, N-M, n, lower.tail=F)

  if (neeat.par$verbose)
    c(z.score, p.value, avg.score, var.score, m)
  else
    c(z.score, p.value)
}

