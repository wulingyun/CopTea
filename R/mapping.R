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

  id1 <- mapping21(symbols, db.SYMBOL)
  symbols <- symbols[!(symbols %in% names(id1))]
  id2 <- mapping21(symbols, db.ALIAS)
  symbols <- symbols[!(symbols %in% names(id2))]
  
  if (length(symbols) > 0)
  {
    warning("The following symbols can not be mapped: ", paste(symbols, collapse=" "))
  }

  c(id1, id2)
}


#' Mapping central id to gene symbols
#'
#'
#' @export
id_symbols <- function(ids, species = "Human")
{
  if (species == "Human") {
    require(org.Hs.eg.db)
    db.SYMBOL <- org.Hs.egSYMBOL
  }
  else if (species == "Mouse") {
    require(org.Mm.eg.db)
    db.SYMBOL <- org.Mm.egSYMBOL
  }
  else {
    stop("Unsupported species!")
  }
  
  symbols <- mapping12(ids, db.SYMBOL)
  ids <- ids[!(ids %in% names(symbols))]
  
  if (length(ids) > 0)
  {
    warning("The following ids can not be mapped: ", paste(ids, collapse=" "))
  }
  
  symbols
}


mapping12 <- function(id, db)
{
  map <- as.matrix(toTable(db))
  match <- map[,1] %in% id
  id <- map[match, 2]
  names(id) <- map[match, 1]
  id
}


mapping21 <- function(id, db)
{
  map <- as.matrix(toTable(db))
  match <- map[,2] %in% id
  id <- map[match, 1]
  names(id) <- map[match, 2]
  id
}
