#' Mapping gene symbols to central id
#'
#'
#' @export
id_mapping <- function(symbols, species = "Human")
{
  if (species == "Human") {
    require(org.Hs.eg.db)
    db.SYMBOL <- org.Hs.egSYMBOL2EG
    db.ALIAS <- org.Hs.egALIAS2EG
  }
  else if (species == "Mouse") {
    require(org.Mm.eg.db)
    db.SYMBOL <- org.Mm.egSYMBOL2EG
    db.ALIAS <- org.Mm.egALIAS2EG
  }
  else {
    stop("Unsupported species!")
  }

  id1 <- mapping(symbols, db.SYMBOL)
  symbols <- symbols[!(symbols %in% names(id1))]
  id2 <- mapping(symbols, db.ALIAS)
  symbols <- symbols[!(symbols %in% names(id2))]
  
  if (length(symbols) > 0)
  {
    warning("The following symbols can not be mapped: ", paste(symbols, collapse=" "))
  }

  c(id1, id2)
}


mapping <- function(id, db)
{
  map <- as.matrix(toTable(db))
  match <- map[,2] %in% id
  id <- map[match, 1]
  names(id) <- map[match, 2]
  id
}
