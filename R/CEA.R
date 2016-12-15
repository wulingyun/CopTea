#' CEA
#' 
#' A novel combination-based method for gene set functional enrichment analysis.
#' 
#' This function performs the CEA algorithm to perform the enrichment analysis.
#' CEA is based on a multi-objective optimization problem, 
#' and the adapted IMPROVED GREEDY algorithm was used to approximatively solve the problem.
#' 
#' @param resultmap The annotation matrix of enrichment analysis. Rows are genes and columns are categories.
#' @param active_gene The input character active gene list.
#' @param d A stochastic tolerance paramter to add randomness to CEA algorithm.
#' @param times The repeated time of CEA algorithm. Positive integer variable.
#' @param trace Logical variable indicated whether tracing information on the solving progress is produced.
#'
#' @return This function will return a list with the following components:
#'   \item{p.values}{The sorted Fisher's exact test p-values of the identified category sets.}
#'   \item{coverage}{The related coverage of the corresponding category set.}
#'   \item{category}{The identified category set corresponding to the sorted p-values.}
#'   \item{annotation}{The pre-processed original annotation matrix.}
#' 
#' @references Duan-Chen Sun, Yin-Liang Liu, Miao Li, Xiang-Sun Zhang, Ling-Yun Wu.
#' CEA: A novel combination-based method for gene set functional enrichment analysis
#' Manuscript, 2016.
#' 
#' @export

CEA <- function(resultmap, active_gene, d = 0, times = 1, trace = TRUE){

	gene <- rownames(resultmap)  						
	UG <- gene %in% active_gene 						# An indicator vector shows the genes need to be covered.
	UO <- !UG  											        # UG and UO is complementary.
	coverage <- colSums(resultmap[UG, , drop=FALSE])
	M_search <- resultmap[, coverage != 0, drop=FALSE] 	# Restrict the candidate categories on the terms annotating at least one genes in UG.
	Ne <- rowSums(M_search) 							  # The annotated number of each gene.
	if(sum(Ne == 0) > 0){     							# There may exist several genes not covered by any categories due to previous preprocessing.
		index <- which(Ne == 0)
		M_search <- M_search[-index, , drop=FALSE]  	# Delete these genes.
		UG <- UG[-index, drop=FALSE]
		UO <- !UG
		Ne <- Ne[-index]
	}
	active_gene <- intersect(rownames(M_search), active_gene)
	m <- nrow(M_search)  								  # The number of genes.
	n <- ncol(M_search)  							  	# The number of categories.
	p <- length(active_gene)							# The number of active genes.
	gene <- rownames(M_search)						# The gene symbols of all covered genes.
	print(paste("Number of annotated active genes:", p))
	
	
	result <-list()									    	# The final result of CEA method.
	result[[1]] <- NULL									  # P values of the final results.
	result[[2]] <- NULL								  	# Coverages of the final results.
	result[[3]] <- list()							  	# Categories of the final results.
	names(result) <- c("p_values", "Coverage", "Category")
	interval <- 0
	for (v in 1:times){
		ind_gene <- rep(0, m)    						# Indicator vector of selected active genes.
		ind_term <- rep(0, n)   						# Indicator vector of selected candidate categories.
		names(ind_gene) <- gene
		names(ind_term) <- colnames(M_search)  

		i <- 1
		M_search1 <- M_search
		weight <- colSums(M_search[UO, ]/ Ne[UO]) + 1e-10					# w(S) vector for each category.
		p_values <- NULL
		Coverage <- NULL
		Category <- list()
		while(sum(ind_gene) < p){
			addnew <- colSums(M_search1[UG & !ind_gene, , drop=FALSE])		# |(intersect(S,UG))\C|
			sigma  <- weight / addnew
			index  <- which(sigma <= (1 + d) * min(sigma))
			index  <- sample(c(index, index), 1)  							# One useful skill.
			ind_gene[names(which(M_search1[UG, index] == 1))] <- 1			# Update the covered active gene and selected category set.
			ind_term[names(index)] <- 1
			Coverage[i]   <- sum(ind_gene & UG) / p
			Category[[i]] <- names(which(ind_term == 1))
			function_gene <- names(which(rowSums(M_search[, Category[[i]], drop=FALSE]) > 0))
			s <- length(function_gene)
			k <- length(intersect(active_gene,function_gene))				# The intersection set of active gene and a GO term.
			p_values[i] <- phyper(k-1,s, m-s, p, lower.tail=FALSE)
			
			cover2 <- colSums(M_search1[UG & !ind_gene, , drop=FALSE])
			index1 <- which(cover2 == 0)
			weight <- weight[-index1]
			M_search1 <- M_search1[,-index1,drop=FALSE]
			i <- i + 1
		}
		if (trace == TRUE){
		  print(paste("Size of categories set in ", v," out of ", times, " repeats : ",i-1, sep=""))
	  }
		interval <- (max(interval) + 1) : (max(interval) + i-1)
		result[[1]][interval] <- p_values
		result[[2]][interval] <- Coverage
		result[[3]][interval] <- Category
	}
	
	print(paste("Number of the identified category sets before unique is",length(result[[1]]),sep=" "))
	tmp <- duplicated(result[[3]])
	if (sum(tmp) > 0){
		result[[1]] <- result[[1]][!tmp]
		result[[2]] <- result[[2]][!tmp]
		result[[3]] <- result[[3]][!tmp]
	}
	print(paste("Number of the identified category sets after  unique is",length(result[[3]]),sep=" "))
	
	
	# Part4: Sort the category sets based on p-values and output the final sorted results.
	ix <- order(result[[1]])
	result[[1]] <- result[[1]][ix]							# P values of the final results.
	result[[2]] <- result[[2]][ix]							# Coverages of the final results.
	result[[3]] <- result[[3]][ix]							# Categories of the final results.

	return(list(p.values = result[[1]], coverage = result[[2]], category = result[[3]], annotation = M_search))
}
