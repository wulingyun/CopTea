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
