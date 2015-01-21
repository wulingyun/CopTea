#' Combination enrichment analysis tool
#' 
#' Functional combination enrichment anslysis of gene sets
#' 
#' This is a novel functional enrichment analysis tool based on combined functional terms.
#' 
#' @param core.sets Logical matrix indicated the core genes associated with specific functions or pathways.
#'  The rows correspond to genes, while the columns represent the functions or pathways.
#' @param gene.set Logical vector indicated the gene set for evaluating. 
#' @param method The method to solve the optimization problem: "MIP" - mixed integer programming (MIP) model 1;
#'  "MIP2" - MIP model 2; "LP" - linear programming relaxation and rounding (LPRR) model 1; "LP2" - LPRR model 2.
#' @param verbose The message output level: 0 - no output; 1 - warning and error only;
#'  2 - normal output; 3 - full output; 4 - debug output. Default: 1;
#'
#' @return This function will return a list with two components:
#'   \item{\code{solution}}{The list of solutions for each coverage level.}
#'   \item{\code{objvalue}}{The list of objective values for each coverage level.}
#'
#' @seealso \code{\link{get_core_sets}}
#' 
#' @import Matrix glpkAPI
#'
#' @export
ceat <- function(core.sets, gene.set, method="LP1", verbose=1)
{
  A <- which(core.sets != 0, arr.ind = T)
  nA <- dim(A)[1]
  nG <- sum(gene.set)
  
  N <- dim(core.sets)[1]
  M <- dim(core.sets)[2]

  solution <- list()
  objvalue <- list()
  
  prob <- initProbGLPK()
  setProbNameGLPK(prob, "CEAT")
  setObjDirGLPK(prob, GLP_MIN)

  if (method == "MIP")
  {
    addRowsGLPK(prob, 2*N+2)
    setRowsNamesGLPK(prob, 1:(2*N+2), c(paste("C1_", 1:N, ""), paste("C2_", 1:N, ""), "C3", "C4"))
    setRowsBndsGLPK(prob, 1:(2*N+2), rep(0, 2*N+2), rep(0, 2*N+2), c(rep(GLP_LO, N), rep(GLP_UP, N), GLP_LO, GLP_UP))
    
    addColsGLPK(prob, M+N)
    setColsNamesGLPK(prob, 1:(M+N), c(paste("X_", 1:M, ""), paste("Y_", 1:N, "")))
    setColsKindGLPK(prob, 1:(M+N), rep(GLP_BV, M+N))
    setColsBndsObjCoefsGLPK(prob, 1:(M+N), rep(0, M+N), rep(1, M+N), c(rep(1/(M^2), M), 1-gene.set))
    
    ia <- c(A[,1], 1:N,         A[,1]+N, (N+1):(N+N), rep(2*N+1, N), rep(2*N+2, N))
    ja <- c(A[,2], (M+1):(M+N), A[,2],   (M+1):(M+N), (M+1):(M+N),   (M+1):(M+N))
    aa <- c(rep(1, nA), rep(-1, N), rep(1, nA), rep(-M, N), gene.set, 1-gene.set)
    loadMatrixGLPK(prob, length(aa), ia, ja, aa)
    
    setMIPParmGLPK(MSG_LEV, verbose)
    setMIPParmGLPK(PRESOLVE, GLP_ON)
    
    B <- N-nG
    for (i in nG:1)
    {
      setRowBndGLPK(prob, 2*N+1, GLP_LO, i, i)
      setRowBndGLPK(prob, 2*N+2, GLP_UP, B, B)
      solveMIPGLPK(prob)
      solution[[i]] <- mipColsValGLPK(prob)
      objvalue[[i]] <- mipObjValGLPK(prob)
      B <- (1-gene.set) %*% solution[[i]][(M+1):(M+N)]
    }
  }
  else if (method == "LP")
  {
    addRowsGLPK(prob, 2*N+1)
    setRowsNamesGLPK(prob, 1:(2*N+1), c(paste("C1_", 1:N, ""), paste("C2_", 1:N, ""), "C3"))
    setRowsBndsGLPK(prob, 1:(2*N+1), rep(0, 2*N+1), rep(0, 2*N+1), c(rep(GLP_LO, N), rep(GLP_UP, N), GLP_LO))
    
    addColsGLPK(prob, M+N)
    setColsNamesGLPK(prob, 1:(M+N), c(paste("X_", 1:M, ""), paste("Y_", 1:N, "")))
    setColsKindGLPK(prob, 1:(M+N), rep(GLP_CV, M+N))
    setColsBndsObjCoefsGLPK(prob, 1:(M+N), rep(0, M+N), rep(1, M+N), c(rep(1/(M^2), M), 1-gene.set))
    
    ia <- c(A[,1], 1:N,         A[,1]+N, (N+1):(N+N), rep(2*N+1, N))
    ja <- c(A[,2], (M+1):(M+N), A[,2],   (M+1):(M+N), (M+1):(M+N))
    aa <- c(rep(1, nA), rep(-1, N), rep(1, nA), rep(-M, N), gene.set)
    loadMatrixGLPK(prob, length(aa), ia, ja, aa)

    setSimplexParmGLPK(MSG_LEV, verbose)

    for (i in nG:1)
    {
      setRowBndGLPK(prob, 2*N+1, GLP_LO, i, i)
      solveSimplexGLPK(prob)
      solution[[i]] <- ceat.rounding(getColsPrimGLPK(prob), core.sets, gene.set, i)
      objvalue[[i]] <- ceat.evaluate(solution[[i]], core.sets, gene.set)
      ceat.check(solution[[i]], core.sets, gene.set, i)
    }
  }
  else if (method == "MIP2")
  {
    addRowsGLPK(prob, N+nA+2)
    setRowsNamesGLPK(prob, 1:(N+nA+2), c(paste("C1_", 1:N, ""), paste("C2_", 1:nA, ""), "C3", "C4"))
    setRowsBndsGLPK(prob, 1:(N+nA+2), rep(0, N+nA+2), rep(0, N+nA+2), c(rep(GLP_LO, N), rep(GLP_UP, nA), GLP_LO, GLP_UP))
    
    addColsGLPK(prob, M+N)
    setColsNamesGLPK(prob, 1:(M+N), c(paste("X_", 1:M, ""), paste("Y_", 1:N, "")))
    setColsKindGLPK(prob, 1:(M+N), rep(GLP_BV, M+N))
    setColsBndsObjCoefsGLPK(prob, 1:(M+N), rep(0, M+N), rep(1, M+N), c(rep(1/(M^2), M), 1-gene.set))
    
    ia <- c(A[,1], 1:N,         (1:nA)+N, (1:nA)+N, rep(N+nA+1, N), rep(N+nA+2, N))
    ja <- c(A[,2], (M+1):(M+N), A[,2],    A[,1]+M,  (M+1):(M+N),    (M+1):(M+N))
    aa <- c(rep(1, nA), rep(-1, N), rep(1, nA), rep(-1, nA), gene.set, 1-gene.set)
    loadMatrixGLPK(prob, length(aa), ia, ja, aa)
    
    setMIPParmGLPK(MSG_LEV, verbose)
    setMIPParmGLPK(PRESOLVE, GLP_ON)
    
    B <- N-nG
    for (i in nG:1)
    {
      setRowBndGLPK(prob, N+nA+1, GLP_LO, i, i)
      setRowBndGLPK(prob, N+nA+2, GLP_UP, B, B)
      solveMIPGLPK(prob)
      solution[[i]] <- mipColsValGLPK(prob)
      objvalue[[i]] <- mipObjValGLPK(prob)
      B <- (1-gene.set) %*% solution[[i]][(M+1):(M+N)]
    }
  }
  else if (method == "LP2")
  {
    addRowsGLPK(prob, N+nA+1)
    setRowsNamesGLPK(prob, 1:(N+nA+1), c(paste("C1_", 1:N, ""), paste("C2_", 1:nA, ""), "C3"))
    setRowsBndsGLPK(prob, 1:(N+nA+1), rep(0, N+nA+1), rep(0, N+nA+1), c(rep(GLP_LO, N), rep(GLP_UP, nA), GLP_LO))
    
    addColsGLPK(prob, M+N)
    setColsNamesGLPK(prob, 1:(M+N), c(paste("X_", 1:M, ""), paste("Y_", 1:N, "")))
    setColsKindGLPK(prob, 1:(M+N), rep(GLP_CV, M+N))
    setColsBndsObjCoefsGLPK(prob, 1:(M+N), rep(0, M+N), rep(1, M+N), c(rep(1/(M^2), M), 1-gene.set))
    
    ia <- c(A[,1], 1:N,         (1:nA)+N, (1:nA)+N, rep(N+nA+1, N))
    ja <- c(A[,2], (M+1):(M+N), A[,2],    A[,1]+M,  (M+1):(M+N))
    aa <- c(rep(1, nA), rep(-1, N), rep(1, nA), rep(-1, nA), gene.set)
    loadMatrixGLPK(prob, length(aa), ia, ja, aa)

    setSimplexParmGLPK(MSG_LEV, verbose)

    for (i in nG:1)
    {
      setRowBndGLPK(prob, N+nA+1, GLP_LO, i, i)
      solveSimplexGLPK(prob)
      solution[[i]] <- ceat.rounding(getColsPrimGLPK(prob), core.sets, gene.set, i)
      objvalue[[i]] <- ceat.evaluate(solution[[i]], core.sets, gene.set)
      ceat.check(solution[[i]], core.sets, gene.set, i)
    }
  }
  
  delProbGLPK(prob)

  list(solution=solution, objvalue=objvalue)
}

ceat.rounding <- function(solution, core.sets, gene.set, alpha)
{
  N <- dim(core.sets)[1]
  M <- dim(core.sets)[2]
  xx <- solution[1:M]
  yy <- solution[(M+1):(M+N)]
  y <- rep(0, length(yy))
  y[order(yy, decreasing=T)[1:alpha]] <- 1
  x <- rep(0, length(xx))
  for (j in which(y == 1))
  {
    x[which.max(core.sets[j, ] * xx)] <- 1
  }
  g_j <- as.numeric(core.sets %*% x)
  y[g_j >= 1] <- 1
  y[g_j < 1] <- 0
  c(x, y)
}

ceat.evaluate <- function(solution, core.sets, gene.set)
{
  N <- dim(core.sets)[1]
  M <- dim(core.sets)[2]
  (1-gene.set) %*% solution[(M+1):(M+N)] + sum(solution[1:M]) / (M^2)
}

ceat.check <- function(solution, core.sets, gene.set, alpha)
{
  solution[solution <= 0.5] <- 0
  solution[solution >= 0.5] <- 1
  N <- dim(core.sets)[1]
  M <- dim(core.sets)[2]
  x <- solution[1:M]
  y <- solution[(M+1):(M+N)]
  yy <- core.sets %*% x
  yy[yy > 1] <- 1
  if (any(y != yy))
  {
    print(paste("Solution is inconsistent at y:", which(y != yy)))
  }
  if (gene.set %*% y < alpha)
  {
    print(paste("Coverage ", alpha, " is not satisfied!"))
  }
}
