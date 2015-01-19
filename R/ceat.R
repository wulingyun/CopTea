#' Combination enrichment analysis tool
#' 
#' Functional combination enrichment anslysis of gene sets
#' 
#' This is a novel functional enrichment analysis tool based on combined functional terms.
#' 
#'
#' @seealso \code{\link{get_core_sets}}
#' 
#' @import Matrix glpkAPI
#'
#' @export
ceat <- function(core.sets, gene.set, method="LP1")
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
    setColsKindGLPK(prob, 1:(M+N), c(rep(GLP_CV, M), rep(GLP_BV, N)))
    setColsBndsObjCoefsGLPK(prob, 1:(M+N), rep(0, M+N), rep(1, M+N), c(rep(1/(M^2), M), 1-gene.set))
    
    ia <- c(A[,1], 1:N,         A[,1]+N, (N+1):(N+N), rep(2*N+1, N), rep(2*N+2, N))
    ja <- c(A[,2], (M+1):(M+N), A[,2],   (M+1):(M+N), (M+1):(M+N),   (M+1):(M+N))
    aa <- c(rep(1, nA), rep(-1, N), rep(1, nA), rep(-M, N), gene.set, 1-gene.set)
    loadMatrixGLPK(prob, length(aa), ia, ja, aa)
    
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
  else if (method == "LP1")
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

    setSimplexParmGLPK(MSG_LEV, GLP_MSG_ERR)

    for (i in nG:1)
    {
      setRowBndGLPK(prob, 2*N+1, GLP_LO, i, i)
      solveSimplexGLPK(prob)
      print(i)
      solution[[i]] <- ceat.rounding(getColsPrimGLPK(prob), core.sets, gene.set, i)
      objvalue[[i]] <- ceat.evaluate(solution[[i]], core.sets, gene.set)
      ceat.check(solution[[i]], core.sets, gene.set, i)
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

    setSimplexParmGLPK(MSG_LEV, GLP_MSG_ERR)

    for (i in nG:1)
    {
      setRowBndGLPK(prob, N+nA+1, GLP_LO, i, i)
      solveSimplexGLPK(prob)
      print(i)
      solution[[i]] <- ceat.rounding(getColsPrimGLPK(prob), core.sets, gene.set, i)
      objvalue[[i]] <- ceat.evaluate(solution[[i]], core.sets, gene.set)
      ceat.check(solution[[i]], core.sets, gene.set, i)
    }
  }

  list(solution=solution, objvalue=objvalue)
}

ceat.rounding <- function(solution, core.sets, gene.set, alpha)
{
  N <- dim(core.sets)[1]
  M <- dim(core.sets)[2]
  xx <- solution[1:M]
  yy <- solution[(M+1):(M+N)]
  a_g <- 1 / (sum(gene.set) - alpha + 1)
  y <- yy
  y[y >= a_g] <- 1
  y[y < a_g] <- 0
  f_j <- yy / rowSums(core.sets)
  f_j[y < 1] <- Inf
  f_i <- colSums(t(t(core.sets) * xx) >= f_j)
  x <- xx
  x[f_i >= 1] <- 1
  x[f_i < 1] <- 0
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
