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
#' @param method A string indicated the NEEAT model, including "neeat", "neeat_h1", "neeat_h2" and "hyper". 
#' Use "hyper" for traditional hypergeometric test, in which the network information is ignored.
#' @param max.depth Integer for the maximum depth considered in the NEEAT models.
#' @param rho The weight parameter for depths.
#' @param n.perm The number of permutations for calculating p-values. Minimum value 100 is required for the model "neeat".
#' Ignored in the models "neeat_h1", "neeat_h2" and "hyper".
#' @param verbose Logical variable indicated whether output the full statistics, such as avg.score, var.score and raw.score.
#' @param n.cpu The number of CPUs/cores used in the parallel computation.
#' @param perm.batch The desired number of permutations in each batch of the parallel computation.
#' 
#' @return This function will return a 2-dimensional array of dimensions \code{c(5, dim(func.gene.sets)[2])},
#' and each column \code{[,j]} containing the following components for the corresponding functional gene set:
#'   \item{\code{z.score}}{The Z-score of the functional gene set}
#'   \item{\code{p.value}}{The statistic significance of the functional gene set}
#'   \item{\code{avg.score}}{The average score in random permutations}
#'   \item{\code{var.score}}{The variance of scores in random permutations}
#'   \item{\code{raw.score}}{The raw score of the functional gene set}
#'
#' @seealso \code{\link{get_func_gene_sets}}, \code{\link{neeat_depths}}
#' 
#' @import Matrix parallel
#'
#' @export
neeat <- function(eval.gene.set, func.gene.sets, net,
                  method = "neeat", max.depth = 1, rho = 1.0, n.perm = 100000,
                  verbose = FALSE, n.cpu = 1, perm.batch = 5000)
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
    scores <- as.numeric(neeat.par$rho[depths + 2] %*% func.gene.sets)
    fgs.sel <- scores > 0
    fgs <- func.gene.sets[, fgs.sel, drop=F]
    raw.score <- scores[fgs.sel]
    if (n.cpu > 1) {
      if (.Platform$OS.type == "windows")
        cl <- makeCluster(n.cpu)
      else
        cl <- makeForkCluster(n.cpu)
      n.cpu <- length(cl)
      neeat.par$n.perm <- ceiling(n.perm / n.cpu)
      perm.score <- clusterCall(cl, neeat_perm, eval.gene.set, fgs, net, raw.score, neeat.par)
      stopCluster(cl)
      perm.score <- matrix(rowMeans(matrix(unlist(perm.score), ncol=n.cpu)), nrow=3)
    }
    else {
      perm.score <- neeat_perm(eval.gene.set, fgs, net, raw.score, neeat.par)
    }
    z.score <- (raw.score - perm.score[2, ]) / sqrt(perm.score[3, ])
    result <- matrix(c(NA, 1.0, NA, NA, 0.0), nrow = 5, ncol = dim(func.gene.sets)[2])
    result[, fgs.sel] <- rbind(z.score, perm.score, raw.score)
  }
  else if (method == "neeat_h1") {
    egs <- Matrix(neeat_depths(eval.gene.set, net, max.depth) >= 0, sparse=T)
    result <- neeat_hyper(egs, func.gene.sets, neeat.par)
  }
  else if (method == "neeat_h2") {
    fgs <- Matrix(neeat_depths(func.gene.sets, net, max.depth) >= 0, sparse=T)
    result <- neeat_hyper(eval.gene.set, fgs, neeat.par)
  }
  else if (method == "hyper") {
    result <- neeat_hyper(eval.gene.set, func.gene.sets, neeat.par)
  }
  else {
    stop("Incorrect method!")
  }
  
  output.names <- c("z.score", "p.value", "avg.score", "var.score", "raw.score")
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


neeat_hyper <- function(eval.gene.set, func.gene.sets, neeat.par)
{
  N <- dim(func.gene.sets)[1]
  M <- colSums(func.gene.sets)
  n <- sum(eval.gene.set)
  jobs <- seq_len(dim(func.gene.sets)[2])
  if (neeat.par$n.cpu > 1) {
    if (.Platform$OS.type == "windows")
      cl <- makeCluster(neeat.par$n.cpu)
    else
      cl <- makeForkCluster(neeat.par$n.cpu)
    result <- parSapply(cl, jobs, neeat_hyper_i, N, n, M, eval.gene.set, func.gene.sets, neeat.par)
    stopCluster(cl)
  }
  else {
    result <- sapply(jobs, neeat_hyper_i, N, n, M, eval.gene.set, func.gene.sets, neeat.par)
  }
  result
}


neeat_hyper_i <- function(i, N, n, M, eval.gene.set, func.gene.sets, neeat.par)
{
  m <- sum(column(func.gene.sets, i) & eval.gene.set)

  if (m > 0) {
    M <- M[i]
    avg.score <- n * M / N
    var.score <- n * (M / N) * ((N - M) / N) * ((N - n) / (N - 1))
    if (var.score > 0)
      z.score <- (m - avg.score) / sqrt(var.score)
    else
      z.score <- NA
    p.value <- phyper(m-1, M, N-M, n, lower.tail=F)
    c(z.score, p.value, avg.score, var.score, m)
  }
  else {
    c(NA, 1.0, NA, NA, 0.0)
  }
}

