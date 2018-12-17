#' Turn a two-columns mapping table to matrix
#'
#' Generate the sparse mapping matrix from a two-columns table.
#'
#' This function generates the sparse logical matrix from the two-columns table returned by \code{\link[AnnotationDbi]{toTable}} 
#' in Bioconductor package \code{AnnotationDbi}.
#' 
#' @param table The two-columns table as matrix or data frame. The other columns will not be used, if available.
#' @param rows The row names for the mapping matrix.
#' @param cols The column names for the mapping matrix.
#' 
#' @return This function returns a sparse logical matrix.
#' 
#' @seealso \code{\link[AnnotationDbi]{toTable}}
#'
#' @examples
#' 
#' \dontrun{
#' if (!requireNamespace("BiocManager"))
#'   install.packages("BiocManager")
#' BiocManager::install("org.Hs.eg.db")
#' x <- toMatrix(toTable(org.Hs.egGO2ALLEGS))
#' }
#'
#' @import Matrix
#' 
#' @export
toMatrix <- function(table, rows = unique(table[,1]), cols = unique(table[,2]))
{
  rid <- seq_along(rows)
  names(rid) <- rows
  cid <- seq_along(cols)
  names(cid) <- cols
  table <- table[table[,1] %in% rows & table[,2] %in% cols, ]
  sparseMatrix(rid[table[,1]], cid[table[,2]], x = T, dims = c(length(rows), length(cols)), dimnames = list(rows, cols))
}


#' Extract a column from a matrix
#' 
#' Extract a specified column from a sparse matrix rapidly
#' 
#' This function use faster extraction algorithm for the \code{\link[=CsparseMatrix-class]{CsparseMatrix}} class in the package \pkg{Matrix}.
#' 
#' @param m The matrix
#' @param i The column index
#' 
#' @return This function will return the specified column as a vector of corresponding type.
#' 
#' @import Matrix
#' 
#' @export
column <- function(m, i)
{
  if (inherits(m, "CsparseMatrix")) {
    v <- vector(typeof(m@x), m@Dim[1])
    p <- (m@p[i]+1):m@p[i+1]
    if (p[1] <= p[length(p)])
      v[m@i[p]+1] <- m@x[p]
  }
  else
    v <- m[,i]
  v
}


#' Extract a sub-matrix from a matrix
#' 
#' Extract a specified sub-matrix from a sparse matrix rapidly
#' 
#' This function use faster extraction algorithm for the \code{\link[=CsparseMatrix-class]{CsparseMatrix}} class in the package \pkg{Matrix}.
#' 
#' @param m The matrix
#' @param rows The row indexes
#' @param cols The column indexes
#' 
#' @return This function will return the specified sub-matrix as a matrix of corresponding type.
#' 
#' @import Matrix
#' 
#' @export
submatrix <- function(m, rows, cols)
{
  sapply(cols, function(i) column(m, i)[rows])
}


nnzero <- function(m, r, c)
{
  in.r <- if (missing(r)) function(i) T else function(i) r[i]
  in.c <- if (missing(c)) function(i) T else function(i) c[i]
  fun <- function(i) if (m@p[i] < m@p[i+1]) (m@p[i]+1):m@p[i+1] else NULL
  if (sum(r) == 0 || sum(c) == 0)
    0
  else if (inherits(m, "CsparseMatrix")) {
    p <- unlist(lapply(which(c), fun))
    sum(m@x[p[in.r(m@i[p]+1)]] != 0)
  }
  else if (inherits(m, "RsparseMatrix")) {
    p <- unlist(lapply(which(r), fun))
    sum(m@x[p[in.c(m@i[p]+1)]] != 0)
  }
  else
    Matrix::nnzero(m[r,c])
}

