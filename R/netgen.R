#' NetGen
#' 
#' A novel enhanced network-based generative model for gene set functional enrichment analysis
#' 
#' This function perform the gene set enrichment analysis using network-based generative model, NetGen. 
#' An additional protein-protein interaction (PPI) network was explicitly used to assist the functional analysis.
#' A greedy-based approximate algorithm was peformed to seek for a sub-optimal solution of the log-likelihood function.
#' 
#' @param annotation The input annotation matrix. Each row is a gene and each column is a functional category.
#' @param PPI The adjacent matrix of biological network.
#' @param active_gene The input character active gene list.
#' @param p1 The probability of the core genes to be activated. 
#' If \code{p1} is a numeric value, \code{netgen} perform the enrichment analysis using the given parameter.
#' If \code{p1} is a numeric vector, \code{netgen} perform the enrichment analysis basing on a mixed parameter selection strategy. 
#' @param p2 The probability of the peripheral genes to be activated.
#' Both numeric value and numeric vector are acceptable as introduced in \code{p1}.
#' @param q The probability of all the other genes to be activated by external factor, such as noise, uncontrollable error in experiment, the incompletion of GO annotation and so on.
#' Both numeric value and numeric vector are acceptable as introduced in \code{p1}.
#' @param alpha A positive number to balance the log-likelihood and penalization term.
#' Pre-set value is 3.
#' @param trace Logical variable indicated whether tracing information on the solving progress is produced.
#' 
#' @return The returns of this function \code{fresult} depend on the model input parameter \code{p1}, \code{p2} and \code{q}. 
#' A matrix of enriched categories and its corresponding Fisher's exact test p-value is returned if the the model input parameter \code{p1}, \code{p2} and \code{q} are all numeric values.
#' A list is returned if at least one of the model input parameter \code{p1}, \code{p2} and \code{q} is numeric vector.
#' Each element of the list is a parameter combination result matrix and a combined p-value.
#' 
#' @references Duanchen Sun, Yinliang Liu, Xiang-Sun Zhang, Ling-Yun Wu.
#' NetGen: a novel network-based generative model for gene set functional enrichment analysis.
#' Manuscript, 2016.
#' 
#' @import Matrix
#' 
#' @export
netgen <- function(annotation, PPI, active_gene, p1, p2, q, alpha = 3, trace=TRUE){
	# NetGen: a novel network-based generative model for gene set functional enrichment analysis.
	
	if (length(p1) == 1 && length(p2) == 1 && length(q) == 1){
		print('Compute the enriched categories using the given parameter combination.')
		Enriched_term <- heuristic(annotation, PPI, active_gene, p1, p2, q, alpha, indicator = trace)
		pvalue <- NULL														# Compute the p-value of identified term using hyper-geometric distribution.
		for (j in 1:length(Enriched_term)){
			function_gene <- rownames(annotation)[which(annotation[,Enriched_term[j]] == 1)]
			N <- nrow(annotation)											# All participant genes.
			M <- length(function_gene)										# Gene number in a same GO term, which owns the same function.
			n <- length(active_gene)										# Active gene number.
			k <- length(intersect(active_gene,function_gene))				# The intersection set of active gene and a GO term.
			pvalue[j] <- phyper(k-1,M, N-M, n, lower.tail=FALSE)
		}
		sp <- sort(pvalue,index.return=T)
		tmp_result <- matrix(0,length(Enriched_term),2,dimnames=list(1:length(Enriched_term),c("GO ID","p-value")))
		tmp_result[,1] <- Enriched_term[sp$ix]
		tmp_result[,2] <- pvalue[sp$ix]
		fresult <- as.data.frame(tmp_result)
	}else{
		print('Compute the enriched categories using a mixed parameter selection strategy.')
		
		sums <- length(p1) * length(p2) * length(q)
		p_pool <- matrix(0,sums,3)												# All candidate parameter combinations.
		v <- 1
		for (i in 1:length(p1)){
			for (j in 1:length(p2)){
				for (k in 1:length(q)){
					p_pool[v,] <- c(p1[i],p2[j],q[k])
					v <- v + 1
				}
			}
		}
		result <- list()
		combined_p <- NULL
		for(i in 1:sums){
			print(paste("Computing parameter combination p1=", p_pool[i,1], ", p2=", p_pool[i,2], ", q=", p_pool[i,3], sep=''))
			tmp_p1 <- p_pool[i,1]
			tmp_p2 <- p_pool[i,2]
			tmp_q  <- p_pool[i,3]
			Enriched_term <- heuristic(annotation, PPI, active_gene, tmp_p1, tmp_p2, tmp_q, alpha, indicator = trace)
			pvalue <- NULL														# Compute the p-value of identified term using hyper-geometric distribution.
			for (j in 1:length(Enriched_term)){
				function_gene <- rownames(annotation)[which(annotation[,Enriched_term[j]] == 1)]
				N <- nrow(annotation)											# All participant genes.
				M <- length(function_gene)										# Gene number in a same GO term, which owns the same function.
				n <- length(active_gene)										# Active gene number.
				k <- length(intersect(active_gene,function_gene))				# The intersection set of active gene and a GO term.
				pvalue[j] <- phyper(k-1,M, N-M, n, lower.tail=FALSE)
			}
			sp <- sort(pvalue,index.return=T)
			tmp_result <- matrix(0,length(Enriched_term),2,dimnames=list(1:length(Enriched_term),c("GO ID","p-value")))
			tmp_result[,1] <- Enriched_term[sp$ix]
			tmp_result[,2] <- pvalue[sp$ix]

			result[[i]] <- as.data.frame(tmp_result)
			names(result)[i] <- paste("p1=",tmp_p1," p2=",tmp_p2," q=",tmp_q,sep="")
			
			function_gene <- rownames(annotation)[rowSums(annotation[,Enriched_term,drop=FALSE]) != 0]
			N <- nrow(annotation)												# All participant genes.
			M <- length(function_gene)											# Gene number in a same GO term, which owns the same function.
			n <- length(active_gene)											# Active gene number.
			k <- length(intersect(active_gene,function_gene))					# The intersection set of active gene and a GO term.
			combined_p[i] <- phyper(k-1,M, N-M, n, lower.tail=FALSE)			# equivalent to phyper(k-1,n, N-n, M, lower.tail=FALSE) 
		}
		fresult <- list()
		fresult[['mix_result']] <- result
		fresult[['Term_combined_pvalue']] <- combined_p
	}
	return(fresult)
}


score <- function(annotation, PPI, active_gene, active_term, p1, p2, q, alpha){
	# Compute the score of the log-likelihood function.
	
	all_gene <- rownames(annotation)
	core_gene <- all_gene[rowSums(annotation[, active_term,drop=FALSE]) != 0]	# Core gene set.
	yes_active_core_gene <- intersect(active_gene, core_gene)					# Active   core gene set.
	non_active_core_gene <- setdiff(core_gene, active_gene)						# Inactive core gene set.
	
	common <- intersect(core_gene, rownames(PPI))								# Note that several core genes are not exist in the PPI network.
	neighbour  <- rownames(PPI)[rowSums(PPI[, common, drop=FALSE]) != 0]		# Search the neighbour node of core gene existing in PPI network. 
	peripheral <- setdiff(neighbour, core_gene)									# Peripheral gene set.

	N1 <- length(yes_active_core_gene)
	E1 <- sum(annotation[non_active_core_gene, active_term])
	N2 <- length(intersect(peripheral, active_gene))
	E2 <- sum(PPI[common, setdiff(peripheral, active_gene)])
	
	colored <- union(peripheral, core_gene)
	remain  <- setdiff(all_gene, colored)
	Oa <- length(intersect(remain, active_gene)) 
	Oi <- length(setdiff(remain, active_gene))

	p <- N1*log(p1) + E1*log(1-p1) + N2*log(p2) + E2*log(1-p2) + Oa*log(q) + Oi*log(1-q) - alpha*length(active_term)
	return(p)
}


heuristic <- function(annotation, PPI, active_gene, p1, p2, q, alpha, indicator=FALSE){
	# Obtain the identified categories using a greedy-based approximate algorithm. The details can be found in GenGO paper.
	
	C0 <- -Inf
	C1 <- 1
	C2 <- 1
	all_term <- colnames(annotation)
	cur_term <- NULL
	while (C0 < C1 || C0 < C2){
		rem_term <- setdiff(all_term, cur_term)
		if (length(cur_term) == 0){
			C0 <- -Inf
		}else{
			C0 <- score(annotation, PPI, active_gene, cur_term, p1, p2, q, alpha)
		}
		tmp_score1 <- NULL											# Compute the delete term.
		if (length(cur_term) > 1){
			for (i in 1:length(cur_term)){
				tmp_term <- cur_term[-i]
				tmp_score1[i] <- score(annotation, PPI, active_gene, tmp_term, p1, p2, q, alpha)
			}
			t1 <- cur_term[which.max(tmp_score1)]
		}else{
			tmp_score1 <- -Inf
		}
		tmp_score2 <- NULL											# Compute the adding term.
		for (i in 1:length(rem_term)){
			tmp_term <- union(rem_term[i], cur_term)
			tmp_score2[i] <- score(annotation, PPI, active_gene, tmp_term, p1, p2, q, alpha)
		}
		t2 <- rem_term[which.max(tmp_score2)]
		C1 <- max(tmp_score1)										# Compute the delete term.				
		C2 <- max(tmp_score2)										# Compute the adding term.	

		if (max(C0, C1, C2) == C1){
			cur_term <- setdiff(cur_term, t1)
		}
		if (max(C0, C1, C2) == C2){
			cur_term <- union(cur_term, t2)
		}
		if (indicator == TRUE){
			print(paste(C0, C1, C2, sep=" "))
		}
	}
	return(cur_term)
}
