RN_iso_select <-
function(Results, gpv_t = 0.01, method = "BH")
{
	if(!is.list(Results) || !all(c("lpv", "gene_status") %in% names(Results)))
	{
		stop("Results must be the output of RN_iso_calc or RNentropy_iso_switch",
			call. = FALSE)
	}

	if(!is.matrix(Results$lpv) || !is.numeric(Results$lpv))
	{
		stop("Results$lpv must be a numeric matrix", call. = FALSE)
	}
	if(is.null(colnames(Results$lpv)) ||
		(nrow(Results$lpv) > 0 && is.null(rownames(Results$lpv))))
	{
		stop("Results$lpv must have row and column names", call. = FALSE)
	}

	if(length(Results$gene_status) != nrow(Results$lpv))
	{
		stop("Results$gene_status must have one value for each row of Results$lpv",
			call. = FALSE)
	}
	if(nrow(Results$lpv) > 0 && !is.null(names(Results$gene_status)) &&
		!identical(names(Results$gene_status), rownames(Results$lpv)))
	{
		stop("names of Results$gene_status must match row names of Results$lpv",
			call. = FALSE)
	}
	if(any(is.na(Results$gene_status)))
	{
		stop("Results$gene_status cannot contain missing values", call. = FALSE)
	}

	if(length(gpv_t) != 1 || !is.numeric(gpv_t) || is.na(gpv_t) ||
		!is.finite(gpv_t) || gpv_t <= 0 || gpv_t > 1)
	{
		stop("gpv_t must be a single p-value greater than 0 and less than or equal to 1",
			call. = FALSE)
	}
	if(length(method) != 1 || !is.character(method) || is.na(method) ||
		!method %in% stats::p.adjust.methods)
	{
		stop("method must be one of p.adjust.methods", call. = FALSE)
	}

	tested <- Results$gene_status == "TESTED"
	eligible <- matrix(tested, nrow = nrow(Results$lpv),
		ncol = ncol(Results$lpv)) & is.finite(Results$lpv)
	if(any(Results$lpv[eligible] < 0))
	{
		stop("finite LPVs for TESTED genes must be non-negative", call. = FALSE)
	}

	gene.names <- rownames(Results$lpv)
	gene.pvalues <- rep(NA_real_, nrow(Results$lpv))
	for(gene in which(rowSums(eligible) > 0))
	{
		gene.pvalues[gene] <- .RN_simes(10 ^ -Results$lpv[gene, eligible[gene, ]])
	}
	family <- !is.na(gene.pvalues)

	Results$gpv <- -log10(gene.pvalues)
	Results$gpv_adj <- rep(NA_real_, length(gene.pvalues))
	if(any(family))
	{
		Results$gpv_adj[family] <- -log10(p.adjust(gene.pvalues[family],
			method = method))
	}
	names(Results$gpv) <- gene.names
	names(Results$gpv_adj) <- gene.names

	selected.rows <- family & Results$gpv_adj >= -log10(gpv_t)

	Results$selected <- data.frame(
		gene_status = unname(Results$gene_status[selected.rows]),
		ISO_GPV = unname(Results$gpv[selected.rows]),
		CORR_ISO_GPV = unname(Results$gpv_adj[selected.rows]),
		Results$lpv[selected.rows, , drop = FALSE],
		row.names = gene.names[selected.rows], check.names = FALSE,
		stringsAsFactors = FALSE)

	Results$selected <- Results$selected[order(-Results$selected$CORR_ISO_GPV,
		-Results$selected$ISO_GPV, seq_len(nrow(Results$selected))), , drop = FALSE]

	return(Results)
}
