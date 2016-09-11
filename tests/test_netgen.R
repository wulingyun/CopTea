library(Matrix)
library(CopTea)
setwd("netgen")
load("annotation.RData")								# The example annotation matrix contains 300 biological processes categories.
load("adjacent_matrix.RData")							# The adjacent matrix of a protein-protein interaction network.
load("active_gene.RData")								# A pre-set 84 generated active gene list derived from the true categories as follows.
options(scipen=0)


# Mode1: Identify the enriched categories using a given parameter combination.
True_categories <- c("GO:0019614", "GO:1903249", "GO:2000506", "GO:0015985", "GO:0071962")
Enriched_Category <- netgen(annotation, adj_matrix, active_gene, p1 = 0.8, p2 = 0.1, q = 0.001, alpha = 5, trace=TRUE)
print('Identified enriched categories via NetGen are')
print(Enriched_Category)
print('False Negative categories identified via NetGen are')
print(setdiff(True_categories, Enriched_Category[,1]))


# Mode2: Identify the enriched categories basing on a mixed parameter selection strategy.
p1 <- c(0.5, 0.8)
p2 <- c(0.1, 0.3)
q  <- 0.001
Enriched_Category <- netgen(annotation, adj_matrix, active_gene, p1, p2, q, alpha = 3, trace=FALSE)
print('The combined p-values of mixed parameter selection strategy are')
print(Enriched_Category$Term_combined_pvalue)
print('The most enriched categories and its corresponding parameter combination are')
print(Enriched_Category$mix_result[which.min(Enriched_Category$Term_combined_pvalue)])
