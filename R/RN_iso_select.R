RN_iso_select <-
function(Results, gene_pv = 0.01)
{
	if(!is.list(Results) || !all(c("pv", "gene_status") %in% names(Results)))
	{
		stop("Results must be the output of RN_iso_calc or RNentropy_iso_switch",
			call. = FALSE)
	}

	if(!is.matrix(Results$pv) || !is.numeric(Results$pv))
	{
		stop("Results$pv must be a numeric matrix", call. = FALSE)
	}
	if(is.null(colnames(Results$pv)) ||
		(nrow(Results$pv) > 0 && is.null(rownames(Results$pv))))
	{
		stop("Results$pv must have row and column names", call. = FALSE)
	}

	if(length(Results$gene_status) != nrow(Results$pv))
	{
		stop("Results$gene_status must have one value for each row of Results$pv",
			call. = FALSE)
	}
	if(nrow(Results$pv) > 0 && !is.null(names(Results$gene_status)) &&
		!identical(names(Results$gene_status), rownames(Results$pv)))
	{
		stop("names of Results$gene_status must match row names of Results$pv",
			call. = FALSE)
	}
	if(any(is.na(Results$gene_status)))
	{
		stop("Results$gene_status cannot contain missing values", call. = FALSE)
	}

	if(length(gene_pv) != 1 || !is.numeric(gene_pv) || is.na(gene_pv) ||
		!is.finite(gene_pv) || gene_pv <= 0 || gene_pv > 1)
	{
		stop("gene_pv must be a single p-value greater than 0 and less than or equal to 1",
			call. = FALSE)
	}

	tested <- Results$gene_status == "TESTED"
	eligible <- matrix(tested, nrow = nrow(Results$pv),
		ncol = ncol(Results$pv)) & is.finite(Results$pv)
	if(any(Results$pv[eligible] < 0))
	{
		stop("finite p-values for TESTED genes must be non-negative", call. = FALSE)
	}

	gene.names <- rownames(Results$pv)
	gene.pvalues <- rep(NA_real_, nrow(Results$pv))
	for(gene in which(rowSums(eligible) > 0))
	{
		gene.pvalues[gene] <- .RN_simes(10 ^ -Results$pv[gene, eligible[gene, ]])
	}
	testable <- !is.na(gene.pvalues)

	Results$gene_pv <- -log10(gene.pvalues)
	names(Results$gene_pv) <- gene.names

	selected.rows <- testable & Results$gene_pv >= -log10(gene_pv)

	Results$selected <- data.frame(
		gene_status = unname(Results$gene_status[selected.rows]),
		ISO_GENE_PV = unname(Results$gene_pv[selected.rows]),
		Results$pv[selected.rows, , drop = FALSE],
		row.names = gene.names[selected.rows], check.names = FALSE,
		stringsAsFactors = FALSE)

	Results$selected <- Results$selected[order(-Results$selected$ISO_GENE_PV,
		seq_len(nrow(Results$selected))), , drop = FALSE]

	return(Results)
}
