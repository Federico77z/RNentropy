RN_iso_select <-
function(Results, lpv_t = 0.01, method = "BH")
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
	if(is.null(rownames(Results$lpv)) || is.null(colnames(Results$lpv)))
	{
		stop("Results$lpv must have row and column names", call. = FALSE)
	}

	if(length(Results$gene_status) != nrow(Results$lpv))
	{
		stop("Results$gene_status must have one value for each row of Results$lpv",
			call. = FALSE)
	}
	if(!is.null(names(Results$gene_status)) &&
		!identical(names(Results$gene_status), rownames(Results$lpv)))
	{
		stop("names of Results$gene_status must match row names of Results$lpv",
			call. = FALSE)
	}
	if(any(is.na(Results$gene_status)))
	{
		stop("Results$gene_status cannot contain missing values", call. = FALSE)
	}

	if(length(lpv_t) != 1 || !is.numeric(lpv_t) || is.na(lpv_t) ||
		!is.finite(lpv_t) || lpv_t <= 0 || lpv_t > 1)
	{
		stop("lpv_t must be a single p-value greater than 0 and less than or equal to 1",
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

	Results$lpv_adj <- matrix(NA_real_, nrow = nrow(Results$lpv),
		ncol = ncol(Results$lpv), dimnames = dimnames(Results$lpv))
	colnames(Results$lpv_adj) <- paste("CORR_ISO_LPV",
		sub("^ISO_LPV_", "", colnames(Results$lpv)), sep = "_")

	if(any(eligible))
	{
		pvalues <- 10 ^ -Results$lpv[eligible]
		Results$lpv_adj[eligible] <- -log10(p.adjust(pvalues, method = method))
	}

	selected.rows <- tested & apply(Results$lpv_adj, 1,
		function(x) any(x >= -log10(lpv_t), na.rm = TRUE))
	selected.names <- rownames(Results$lpv)[selected.rows]

	Results$selected <- data.frame(
		gene_status = unname(Results$gene_status[selected.rows]),
		Results$lpv[selected.rows, , drop = FALSE],
		Results$lpv_adj[selected.rows, , drop = FALSE],
		row.names = selected.names, check.names = FALSE,
		stringsAsFactors = FALSE)

	if(any(selected.rows))
	{
		strongest <- apply(Results$lpv_adj[selected.rows, , drop = FALSE], 1,
			max, na.rm = TRUE)
		Results$selected <- Results$selected[order(-strongest,
			seq_along(strongest)), , drop = FALSE]
	}

	return(Results)
}
