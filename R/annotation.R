#' Get annotations from databases
#'
#' @export
get_annotation <- function(species, STRING.version = "9_1", STRING.threshold = 900)
{
  if (species == "Human") {
    require(org.Hs.eg.db)
    pre <- "hsa"
    taxonomy <- 9606
    db.GO <- org.Hs.egGO2ALLEGS
    db.PROT <- org.Hs.egENSEMBLPROT
  }
  else if (species == "Mouse") {
    require(org.Mm.eg.db)
    pre <- "mmu"
    taxonomy <- 10090
    db.GO <- org.Mm.egGO2ALLEGS
    db.PROT <- org.Mm.egENSEMBLPROT
  }
  else if (species == "Yeast") {
    require(org.Sc.sgd.db)
    pre <- "sce"
    taxonomy <- 4932
    db.GO <- org.Sc.sgdGO2ALLORFS
    db.PROT <- org.Sc.sgdENSEMBLPROT
  }
  else {
    stop("Unsupported species!")
  }
  
  data <- list()
  
  require(KEGG.db)
  map.GO <- as.matrix(toTable(db.GO))
  map.GO <- map.GO[, c(1,2)]
  map.KEGG <- as.matrix(toTable(KEGGPATHID2EXTID))
  map.KEGG <- cbind(map.KEGG[,2], map.KEGG[,1])
  map.KEGG <- map.KEGG[grepl(pre, map.KEGG[,2]), ]
  map <- rbind(map.GO, map.KEGG)

  data$annotations <- toMatrix(map)
  data$genes <- rownames(data$annotations)
  data$terms <- colnames(data$annotations)
  
  map <- as.matrix(toTable(db.PROT))
  
  data$proteins <- unique(map[,2])
  data$gene2protein <- toMatrix(map, data$genes, data$proteins)
  
  require(STRINGdb)
  string.id <- paste(taxonomy, data$proteins, sep=".")
  string.db <- STRINGdb$new(version=STRING.version, species=taxonomy, score_threshold=STRING.threshold)
  ppi <- string.db$get_interactions(string.id)
  ppi <- matrix(c(ppi[,1], ppi[,2], ppi[,2], ppi[,1]), ncol=2)
  ppi <- toMatrix(ppi, string.id, string.id)
  
  data$network <- (data$gene2protein %*% ppi %*% t(data$gene2protein)) > 0
  
  data
}
