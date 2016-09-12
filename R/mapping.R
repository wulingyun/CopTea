#' Mapping gene symbols to central ids
#'
#' Map gene symbols to gene central ids by using organism annotation package in Bioconductor.
#' 
#' @param symbols The vector of gene symbols.
#' @param species The reference species name.
#' @param warning.unmapped Logical variable indicated whether give warnings when 
#' there are gene symbols that can not be mapped.
#' 
#' @return This function will return a vector of gene ids mapped from given gene symbols.
#' 
#' @seealso \code{\link{id_mapping_to_symbol}}, \code{\link{id_mapping}},
#' \code{\link{id_mapping_species}}
#'
#' @export
id_mapping_from_symbol <- function(symbols, species = "Human", warning.unmapped = TRUE)
{
  if (species == "Human") {
    require_bioc("org.Hs.eg.db")
    db.SYMBOL <- org.Hs.eg.db::org.Hs.egSYMBOL2EG
    db.ALIAS <- org.Hs.eg.db::org.Hs.egALIAS2EG
  }
  else if (species == "Mouse") {
    require_bioc("org.Mm.eg.db")
    db.SYMBOL <- org.Mm.eg.db::org.Mm.egSYMBOL2EG
    db.ALIAS <- org.Mm.eg.db::org.Mm.egALIAS2EG
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


#' Mapping central ids to gene symbols
#'
#' Map gene central ids to gene symbols by using organism annotation package in Bioconductor.
#' 
#' @param ids The vector of gene central ids.
#' @param species The reference species name.
#' @param warning.unmapped Logical variable indicated whether give warnings when 
#' there are gene central ids that can not be mapped.
#' 
#' @return This function will return a vector of gene symbols mapped from given gene central ids.
#' 
#' @seealso \code{\link{id_mapping_from_symbol}}, \code{\link{id_mapping}}, 
#' \code{\link{id_mapping_species}}
#'
#' @export
id_mapping_to_symbol <- function(ids, species = "Human", warning.unmapped = TRUE)
{
  if (species == "Human") {
    require_bioc("org.Hs.eg.db")
    db.SYMBOL <- org.Hs.eg.db::org.Hs.egSYMBOL
  }
  else if (species == "Mouse") {
    require_bioc("org.Mm.eg.db")
    db.SYMBOL <- org.Mm.eg.db::org.Mm.egSYMBOL
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


#' Mapping between different gene ids
#'
#' Map gene ids to other gene ids by using organism annotation package in Bioconductor.
#' 
#' @param ids The vector of input gene ids.
#' @param db The gene id mapping database.
#' @param from The index of input gene id in the mapping database.
#' @param to The index of output gene id in the mapping database.
#' 
#' @return This function will return a vector of output gene ids mapped from given gene ids.
#'
#' @seealso \code{\link{id_mapping_from_symbol}}, \code{\link{id_mapping_to_symbol}}, 
#' \code{\link{id_mapping_species}}
#' 
#' @export
id_mapping <- function(ids, db, from = 1, to = 2)
{
  map <- as.matrix(AnnotationDbi::toTable(db))
  ids <- toupper(ids)
  map[, from] <- toupper(map[, from])
  match <- map[, from] %in% ids
  ids.mapped <- map[match, to]
  names(ids.mapped) <- map[match, from]
  ids.unmapped <- ids[!(ids %in% names(ids.mapped))]
  list(mapped=ids.mapped, unmapped=ids.unmapped)
}


#' Mapping gene central ids from one species to other species
#'
#' Map gene central ids in one species to other speicies by using organism annotation package in Bioconductor.
#' 
#' @param ids The vector of input gene central ids.
#' @param species.from The species of input gene ids.
#' @param species.to The species of output gene ids.
#' @param rm.na Logical variable indicated whether remove NA in the output.
#' 
#' @return This function will return a vector of gene central ids mapped from given gene central ids.
#' 
#' @seealso \code{\link{id_mapping_from_symbol}}, \code{\link{id_mapping_to_symbol}}, 
#' \code{\link{id_mapping}}
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
