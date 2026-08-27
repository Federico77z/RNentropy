RNentropy_iso_switch <-
function(file, tr.col, gene.col, design = NULL, header = TRUE, skip = 0,
	skip.col = NULL, col.names = NULL, min_expr = 1, pseudocount = 0.01)
{
	TABLE <- read.table(file, header = header, skip = skip,
		blank.lines.skip = TRUE, comment.char = "#")
	tr.index <- .RN_column_index(TABLE, tr.col, "tr.col")
	gene.index <- .RN_column_index(TABLE, gene.col, "gene.col")

	if(tr.index == gene.index)
	{
		stop("tr.col and gene.col must identify different columns", call. = FALSE)
	}

	skip.index <- integer()
	if(!is.null(skip.col))
	{
		skip.index <- .RN_column_indices(TABLE, skip.col, "skip.col")
		if(tr.index %in% skip.index)
		{
			stop("tr.col cannot be included in skip.col", call. = FALSE)
		}
		if(gene.index %in% skip.index)
		{
			stop("gene.col cannot be included in skip.col", call. = FALSE)
		}
	}

	transcripts <- as.character(TABLE[[tr.index]])
	if(any(is.na(transcripts)) || any(transcripts == "") ||
		anyDuplicated(transcripts))
	{
		stop("transcript identifiers must be univocal and cannot be missing or empty",
			call. = FALSE)
	}

	keep <- setdiff(seq_len(ncol(TABLE)), c(tr.index, skip.index))
	gene.col <- match(gene.index, keep)
	TABLE <- TABLE[keep]
	row.names(TABLE) <- transcripts

	if(!is.null(col.names))
	{
		if(length(col.names) != ncol(TABLE))
		{
			stop("the number of col.names must match the imported columns", call. = FALSE)
		}
		colnames(TABLE) <- col.names
	}

	return(RN_iso_calc(TABLE, gene.col, design, min_expr, pseudocount))
}
