#' Get annotations from databases
#' 
#' Get the annotations needed for functional enrichment analysis from the databases, including GO
#' and KEGG/Reactome pathway annotations, gene to protein mapping, protein interaction network.
#' 
#' @param species A string indicated species name, currently support \code{"Human"}, 
#' \code{"Mouse"} and \code{"Yeast"}.
#' @param filters A vector of string indicated the annotations used, currently support
#' \code{"GO"}, \code{"KEGG"}, \code{"Reactome"}, and \code{"OMIM"} (only valid for
#' Human). 
#' @param STRING.version A string indicated the version number of STRING database used for 
#' extracting protein interaction network.
#' @param STRING.threshold A integer number (0-1000) indicated the threshold for extracting 
#' network from STRING database.
#' 
#' @return This function will return a list with the following components:
#'   \item{anno.matrix}{A logical matrix indicated the annotations, with each row represents 
#'   a gene and each column denotes a term. The row names and column names are set as the 
#'   corresponding gene and term ID respectively.} 
#'   \item{network}{A logical matrix indicated the association (interaction) between genes, 
#'   which is aggregated from the protein association (interaction) network provided by the
#'   STRING database.}
#'   \item{gene2protein}{A logical matrix indicated the mapping from gene to protein, with 
#'   each row represents a gene and each column denotes a protein. The row names and column 
#'   names are set as the corresponding gene and protein ID respectively.}
#'   \item{genes}{A vector of gene IDs.}
#'   \item{terms}{A vector of term IDs.}
#'   \item{proteins}{A vector of protein IDs.}
#'   \item{term.info}{A matrix of term information, with three columns: Term ID, Term category,
#'   and Term description. The row names are set as Term ID.}
#'
#' @export
get_annotations <- function(species, filters = c("GO", "KEGG", "Reactome", "OMIM"),
                            STRING.version = "10", STRING.threshold = 900)
{
  require_bioc("AnnotationDbi")
  if (species == "Human") {
    require_bioc("org.Hs.eg.db")
    species.name <- "Homo sapiens"
    species.pre <- "hsa"
    species.id <- 9606
    db.GO <- org.Hs.eg.db::org.Hs.egGO2ALLEGS
    db.PROT <- org.Hs.eg.db::org.Hs.egENSEMBLPROT
  }
  else if (species == "Mouse") {
    require_bioc("org.Mm.eg.db")
    species.name <- "Mus musculus"
    species.pre <- "mmu"
    species.id <- 10090
    db.GO <- org.Mm.eg.db::org.Mm.egGO2ALLEGS
    db.PROT <- org.Mm.eg.db::org.Mm.egENSEMBLPROT
  }
  else if (species == "Yeast") {
    require_bioc("org.Sc.sgd.db")
    species.name <- "Saccharomyces cerevisiae"
    species.pre <- "sce"
    species.id <- 4932
    db.GO <- org.Sc.sgd.db::org.Sc.sgdGO2ALLORFS
    db.PROT <- org.Sc.sgd.db::org.Sc.sgdENSEMBLPROT
  }
  else {
    stop("Unsupported species!")
  }
  
  data <- list()

  map <- NULL

  if ("GO" %in% filters) {
    if (species != "Human") {
      require_bioc("org.Hs.eg.db")
      map.Human <- as.matrix(AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egGO2ALLEGS))
      map.Human <- map.Human[, c(1,2)]
      genes <- id_mapping_species(map.Human[, 1], "Human", species, FALSE)
      map.Human[, 1] <- genes[map.Human[, 1]]
      map.Human <- map.Human[!is.na(map.Human[, 1]), ]
      map <- rbind(map, map.Human)
    }
    
    map.GO <- as.matrix(AnnotationDbi::toTable(db.GO))
    map.GO <- map.GO[, c(1,2)]
    map <- rbind(map, map.GO)
  }
  
  if ("KEGG" %in% filters) {
    require_bioc("KEGG.db")
    map.KEGG <- as.matrix(AnnotationDbi::toTable(KEGG.db::KEGGPATHID2EXTID))
    map.KEGG <- map.KEGG[, c(2,1)]
    map.KEGG <- map.KEGG[grepl(species.pre, map.KEGG[,2]), ]
    map <- rbind(map, map.KEGG)
  }
  
  if ("Reactome" %in% filters) {
    require_bioc("reactome.db")
    reactome.id <- as.matrix(AnnotationDbi::toTable(reactome.db::reactomePATHID2NAME))
    reactome.id <- reactome.id[grepl(paste("^", species.name, sep=""), reactome.id[, 2]), 1]
    map.REACTOME <- as.matrix(AnnotationDbi::toTable(reactome.db::reactomePATHID2EXTID))
    map.REACTOME <- map.REACTOME[, c(2,1)]
    map.REACTOME <- map.REACTOME[map.REACTOME[, 2] %in% reactome.id, ]
    map.REACTOME[, 2] <- paste("Reactome", map.REACTOME[, 2], sep=":")
    map <- rbind(map, map.REACTOME)
  }

  if ("OMIM" %in% filters && species == "Human") {
    map.OMIM <- as.matrix(AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egOMIM2EG))
    map.OMIM[, 2] <- paste("OMIM", map.OMIM[, 2], sep=":")
    map <- rbind(map, map.OMIM)
  }
  
  data$anno.matrix <- toMatrix(map)
  data$genes <- rownames(data$anno.matrix)
  data$terms <- colnames(data$anno.matrix)
  
  map <- as.matrix(AnnotationDbi::toTable(db.PROT))
  
  data$proteins <- unique(map[,2])
  data$gene2protein <- toMatrix(map, data$genes, data$proteins)
  
  require_bioc("STRINGdb")
  string.id <- paste(species.id, data$proteins, sep=".")
  string.db <- STRINGdb::STRINGdb$new(version=STRING.version, species=species.id, score_threshold=STRING.threshold)
  ppi <- string.db$get_interactions(string.id)
  ppi <- matrix(c(ppi[,1], ppi[,2], ppi[,2], ppi[,1]), ncol=2)
  ppi <- toMatrix(ppi, string.id, string.id)
  
  data$network <- (data$gene2protein %*% ppi %*% t(data$gene2protein)) > 0

  term.info <- NULL
  
  if ("GO" %in% filters) {
    require_bioc("GO.db")
    goterm <- AnnotationDbi::as.list(GO.db::GOTERM)
    goterm <- t(sapply(1:length(goterm), function(i) c(goterm[[i]]@GOID, goterm[[i]]@Ontology, goterm[[i]]@Term)))
    term.info <- rbind(term.info, goterm)
  }
  
  if ("KEGG" %in% filters) {
    require_bioc("KEGG.db")
    kegg <- as.matrix(AnnotationDbi::toTable(KEGG.db::KEGGPATHID2NAME))
    kegg[, 1] <- paste(species.pre, kegg[, 1], sep="")
    kegg <- cbind(kegg[, 1], "KEGG", kegg[, 2])
    term.info <- rbind(term.info, kegg)
  }

  if ("Reactome" %in% filters) {
    require_bioc("reactome.db")
    reactome <- as.matrix(AnnotationDbi::toTable(reactome.db::reactomePATHID2NAME))
    reactome <- reactome[grepl(paste("^", species.name, sep=""), reactome[, 2]), ]
    reactome[, 1] <- paste("Reactome", reactome[, 1], sep=":")
    reactome <- cbind(reactome[, 1], "Reactome", reactome[, 2])
    term.info <- rbind(term.info, reactome)
  }
  
  if ("OMIM" %in% filters && species == "Human") {
    omim <- as.matrix(AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egOMIM2EG))
    omim <- unique(paste("OMIM", omim[, 2], sep=":"))
    omim <- cbind(omim, "OMIM", omim)
    term.info <- rbind(term.info, omim)
  }
  
  term.info <- cbind(term.info, colSums(data$anno.matrix)[term.info[, 1]])
  rownames(term.info) <- term.info[, 1]
  colnames(term.info) <- c("ID", "Category", "Term", "Size")
  data$term.info <- data.frame(term.info, stringsAsFactors=FALSE)
  
  data
}


#' Get significant term list
#'
#' Extract and rank the significant terms with p-values smaller than the given threshold.
#'
#' @param p.values The vector of p-values.
#' @param term.info The matrix of term information with three columns: Term ID, Term category,
#'  and Term description. The row names should be set as Term ID. 
#'  See \code{\link{get_annotations}} for more information.
#' @param threshold The threshold of significant p-value.
#' @param filters The vector of term categories for filtering terms. Set NULL to disable filter.
#' @param z.score The vector of Z-scores associated with p-values, used to sort the terms of same p-values.
#' @param adjust.p The method for adjusting p-values for multiple testing. Use "none" for bypassing.
#' See \code{\link{p.adjust}} for available methods.
#'
#' @return This function returns a matrix with four columns: Rank, Term ID, Term category,
#'  and Term description.
#'  
#' @seealso \code{\link{get_annotations}}
#'
#' @export
get_significant_terms <- function(p.values, term.info, threshold = 0.05, filters = NULL, z.score = NULL, adjust.p = "bonferroni")
{
  p <- data.frame(names(p.values), p.values, stats::p.adjust(p.values, method = adjust.p), stringsAsFactors = F)
  select.term <- p[, 3] <= threshold
  p <- p[select.term, ]
  if (is.null(z.score))
    p <- p[order(p[, 2]), ]
  else
    p <- p[order(p[, 2], -z.score[select.term]), ]
  rank <- cbind(seq_len(nrow(p)), p, term.info[p[, 1], c(2, 3, 4)])
  rownames(rank) <- rank[, 2]
  colnames(rank) <- c("Rank", "ID", "p", adjust.p, "Category", "Term", "Size")
  if (!is.null(filters)) rank <- rank[rank$Category %in% filters, ]
  rank
}


#' Get annotated gene list
#'
#' Extract the genes associated with terms from the annotation matrix.
#'
#' @param anno.matrix A logical matrix indicated the annotations, with each row represents a gene and each column denotes a term.
#'  The row names and column names are set as the corresponding gene and term ID respectively.
#'  See \code{\link{get_annotations}} for more information.
#' @param species A string indicated the species name used for translating gene ID to symbols.
#' @param gene.symbol A logical variable indicated whether use the gene symbol in the output instead of the central ID.
#' 
#' @return This function will return a vector of strings, which is the concatenated IDs (or symbols) of the genes annotated
#' by the corresponding term.
#'  
#' @seealso \code{\link{get_annotations}}
#'
#' @export
get_annotated_genes <- function(anno.matrix, species = "Human", gene.symbol = TRUE)
{
  genes <- rownames(anno.matrix)
  if (gene.symbol) {
    ids <- genes
    symbols <- id_mapping_to_symbol(ids, species)
    genes <- symbols[ids]
    genes[is.na(genes)] <- ids[is.na(genes)]
  }
  terms <- colnames(anno.matrix)
  gene.list <- sapply(seq_along(terms), function(i) paste(genes[anno.matrix[, i]], collapse=", "))
  names(gene.list) <- terms
  gene.list
}


#' Extract GO temrs information from GO annotation
#'
#' Generate the annotation matrix from an object of class "Go3AnnDbBimap" in Bioconductor annotation package.
#'
#' This function generates the annotation matrix required by enrichment analysis methods.
#' 
#' @param go.map An object of class "Go3AnnDbBimap".
#' @param evidence Vector of string to filter the GO annotations, could be "ALL" or a set of evidence codes.
#' @param category Vector of string to filter the GO categories, could be "ALL" or a set of GO categories.
#' @param gene.set Vector of string to filter the genes.
#' 
#' @return This function returns a sparse matrix as required by enrichment analysis methods,
#' with each row represents a gene and each column denotes a term.
#' The row names and column names are set as the corresponding gene and term ID respectively.
#' 
#' @seealso \code{\link{netgen}}
#'
#' @examples
#' 
#' \dontrun{
#' if (!requireNamespace("BiocManager"))
#'   install.packages("BiocManager")
#' BiocManager::install("org.Hs.eg.db")
#' x <- get_GO_annotations(org.Hs.egGO2ALLEGS)
#' }
#'
#' @import Matrix
#' 
#' @export
get_GO_annotations <- function(go.map, evidence = "ALL", category = "ALL", gene.set = NULL)
{
  require_bioc("AnnotationDbi")
  term.table <- AnnotationDbi::toTable(go.map)
  term.table[,1] <- toupper(term.table[,1])
  
  if (evidence != "ALL")
    term.table <- term.table[term.table[,3] %in% evidence, ]
  
  if (category != "ALL")
    term.table <- term.table[term.table[,4] %in% category, ]
  
  if (!is.null(gene.set)) {
    all.gene <- gene.set    
  }
  else
    all.gene <- unique(term.table[,1])
  
  toMatrix(term.table, all.gene)
}


#' Get GO leaf terms
#'
#' Extract the leaf terms from the specific GO domains.
#'
#' @param domains A vector of strings indicated the specific GO domains, possible values include
#'  \code{"BP"}, \code{"CC"}, and \code{"MF"}.
#' 
#' @return This function will return a vector of strings, which is the GO IDs of extracted terms.
#'  
#' @export
get_GO_leafs <- function(domains = c("BP", "CC", "MF"))
{
  require_bioc("AnnotationDbi")
  require_bioc("GO.db")
  terms <- NULL
  if ("BP" %in% domains) {
    x <- is.na(AnnotationDbi::as.list(GO.db::GOBPOFFSPRING))
    terms <- c(terms, names(x)[x])
  }
  if ("CC" %in% domains) {
    x <- is.na(AnnotationDbi::as.list(GO.db::GOCCOFFSPRING))
    terms <- c(terms, names(x)[x])
  }
  if ("MF" %in% domains) {
    x <- is.na(AnnotationDbi::as.list(GO.db::GOMFOFFSPRING))
    terms <- c(terms, names(x)[x])
  }
  terms
}


#' Load a Bioconductor package
#'
#' Load a Bioconductor package, automatically install if need.
#'
#' @param pkg.name The name of package.
#' 
#' @export
require_bioc <- function(pkg.name)
{
  if (!requireNamespace(pkg.name, character.only=T, quietly = T))
  {
    if (!requireNamespace("BiocManager"))
      utils::install.packages("BiocManager")
    BiocManager::install(pkg.name)
    requireNamespace(pkg.name, character.only=T)
  }
}
