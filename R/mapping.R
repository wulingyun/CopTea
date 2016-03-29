#' Mapping gene symbols to central id
#'
#'
#' @export
id_mapping_from_symbol <- function(symbols, species = "Human", warning.unmapped = TRUE)
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

  id1 <- id_mapping(symbols, db.SYMBOL, 2, 1)
  id2 <- id_mapping(id1$unmapped, db.ALIAS, 2, 1)
  
  if (warning.unmapped && length(id2$unmapped) > 0)
  {
    warning("The following symbols can not be mapped: ", paste(id2$unmapped, collapse=" "))
  }

  c(id1$mapped, id2$mapped)
}


#' Mapping central id to gene symbols
#'
#'
#' @export
id_mapping_to_symbol <- function(ids, species = "Human", warning.unmapped = TRUE)
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
  
  id1 <- id_mapping(ids, db.SYMBOL, 1, 2)
  
  if (warning.unmapped && length(id1$unmapped) > 0)
  {
    warning("The following ids can not be mapped: ", paste(id1$unmapped, collapse=" "))
  }
  
  id1$mapped
}


#' Mapping central id to other id
#'
#'
#' @export
id_mapping <- function(id, db, from = 1, to = 2)
{
  map <- as.matrix(toTable(db))
  id <- toupper(id)
  map[, from] <- toupper(map[, from])
  match <- map[, from] %in% id
  id.mapped <- map[match, to]
  names(id.mapped) <- map[match, from]
  id.unmapped <- id[!(id %in% names(id.mapped))]
  list(mapped=id.mapped, unmapped=id.unmapped)
}


#' Mapping gene central id from one species to other species
#'
#'
#' @export
id_mapping_species <- function(ids, species.from, species.to, rm.na = FALSE)
{
  ids <- unique(ids)
  symbols <- id_mapping_to_symbol(ids, species.from, FALSE)
  genes <- id_mapping_from_symbol(symbols, species.to, FALSE)
  genes <- genes[symbols[ids]]
  names(genes) <- ids
  if (rm.na) genes <- genes[!is.na(genes)]
  genes
}
