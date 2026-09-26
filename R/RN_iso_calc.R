RN_iso_calc <-
function(X, gene.col, design = NULL, min_expr = 1, pseudocount = 0.01)
{
	if(!is.data.frame(X))
	{
		X <- as.data.frame(X, stringsAsFactors = FALSE)
	}

	gene.col <- .RN_column_index(X, gene.col, "gene.col")
	rnums <- sapply(X, is.numeric)
	rnums[gene.col] <- FALSE

	if(!any(rnums))
	{
		stop("no numeric expression columns found in X", call. = FALSE)
	}

	genes <- as.character(X[[gene.col]])
	if(any(is.na(genes)) || any(genes == ""))
	{
		stop("gene identifiers cannot be missing or empty", call. = FALSE)
	}

	expression <- X[rnums]
	if(any(is.na(expression)) || any(!is.finite(as.matrix(expression))))
	{
		stop("expression values must be finite and cannot be missing", call. = FALSE)
	}
	if(any(as.matrix(expression) < 0))
	{
		stop("expression values cannot be negative", call. = FALSE)
	}

	.RN_positive_scalar(min_expr, "min_expr", allow.zero = TRUE)
	.RN_positive_scalar(pseudocount, "pseudocount", allow.zero = FALSE)

	if(is.null(design))
	{
		design <- .RN_default_design(ncol(expression))
	}
	.RN_design_check(X, design, rnums)
	replicates <- .RN_get_replicate_list(design)

	gene.names <- sort(unique(genes))
	sample.names <- paste("ISO_PV", colnames(expression), sep = "_")
	if(is.null(colnames(expression)))
	{
		sample.names <- paste("ISO_PV", seq_len(ncol(expression)), sep = "_")
	}

	pv <- matrix(NA_real_, nrow = length(gene.names), ncol = ncol(expression),
		dimnames = list(gene.names, sample.names))
	sample.status <- matrix(NA_character_, nrow = length(gene.names),
		ncol = ncol(expression), dimnames = list(gene.names, sample.names))
	gene.status <- rep("TESTED", length(gene.names))
	names(gene.status) <- gene.names
	gene.rows <- split(seq_len(nrow(X)), factor(genes, levels = gene.names))

	for(gene in gene.names)
	{
		isoforms <- as.matrix(expression[gene.rows[[gene]], , drop = FALSE])

		if(nrow(isoforms) < 2)
		{
			gene.status[gene] <- "NOTEST(ISO)"
			next
		}

		expressed <- isoforms > min_expr
		if(sum(rowSums(expressed) > 0) < 2 || sum(colSums(expressed) > 0) < 2)
		{
			gene.status[gene] <- "NOTEST(EXP)"
			next
		}

		for(sample in seq_len(ncol(isoforms)))
		{
			current <- isoforms[, sample]
			current.total <- sum(current)
			reference.columns <- setdiff(seq_len(ncol(isoforms)), replicates[[sample]])
			reference <- rowSums(isoforms[, reference.columns, drop = FALSE])
			reference.total <- sum(reference)

			if(current.total == 0)
			{
				sample.status[gene, sample] <- "NO_CURRENT_EXPRESSION"
				next
			}
			if(reference.total == 0)
			{
				sample.status[gene, sample] <- "NO_REFERENCE_EXPRESSION"
				next
			}

			current.frequency <- (current / current.total) + pseudocount
			reference.frequency <- (reference / reference.total) + pseudocount
			current.frequency <- current.frequency / sum(current.frequency)
			reference.frequency <- reference.frequency / sum(reference.frequency)

			statistic <- 2 * sum(current * log(current.frequency / reference.frequency))
			pvalue <- pchisq(statistic, nrow(isoforms) - 1, lower.tail = FALSE)
			pv[gene, sample] <- if(pvalue > 0) min(-log10(pvalue), 300) else 300
			sample.status[gene, sample] <- "TESTED"
		}
	}

	res <- data.frame(gene_status = unname(gene.status), pv,
		row.names = gene.names, check.names = FALSE, stringsAsFactors = FALSE)

	return(list(expr = X, design = design, pv = pv,
		gene_status = gene.status, sample_status = sample.status, res = res))
}
