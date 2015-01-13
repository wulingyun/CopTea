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
ceat <- function(core.sets, gene.set)
{
  A <- which(core.sets != 0, arr.ind = T)
  nA <- dim(A)[1]
  
  N <- dim(core.sets)[1]
  M <- dim(core.sets)[2]
  
  prob <- initProbGLPK()
  setProbNameGLPK(prob, "CEAT")
  setObjDirGLPK(prob, GLP_MIN)
  
  addRowsGLPK(prob, 2*N+1)
  addColsGLPK(prob, M+N)
  setRowsNamesGLPK(prob, 1:(2*N+1), c(paste("C1_", 1:N, ""), paste("C2_", 1:N, ""), "C3"))
  setColsNamesGLPK(prob, 1:(M+N), c(paste("X_", 1:M, ""), paste("Y_", 1:N, "")))
  setColsKindGLPK(prob, 1:(M+N), c(rep(GLP_CV, M), rep(GLP_BV, N)))

  setRowsBndsGLPK(prob, 1:(2*N+1), rep(0, 2*N+1), rep(0, 2*N+1), c(rep(GLP_LO, N), rep(GLP_UP, N), GLP_LO))
  setColsBndsObjCoefsGLPK(prob, 1:(M+N), rep(0, M+N), rep(1, M+N), c(rep(1/(M*10), M), 1-gene.set))
  
  ia <- c(A[,1], 1:N,         A[,1]+N, (N+1):(N+N), rep(2*N+1, N))
  ja <- c(A[,2], (M+1):(M+N), A[,2],   (M+1):(M+N), (M+1):(M+N))
  aa <- c(rep(1, nA), rep(-1, N), rep(1, nA), rep(-M, N), gene.set)
  loadMatrixGLPK(prob, length(aa), ia, ja, aa)

  setMIPParmGLPK(PRESOLVE, GLP_ON)

  solution <- list()
  objvalue <- list()
  
  nG <- sum(gene.set)
  for (i in nG:1)
  {
    setRowBndGLPK(prob, 2*N+1, GLP_LO, i, i)
    solveMIPGLPK(prob)
    solution[[i]] <- mipColsValGLPK(prob)
    objvalue[[i]] <- mipObjValGLPK(prob)
  }
  list(solution=solution, objvalue=objvalue)
}
