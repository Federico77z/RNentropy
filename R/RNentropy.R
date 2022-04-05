RNentropy <-
function(file, tr.col, design = NULL, header = TRUE, skip = 0, skip.col = NULL, col.names = NULL)
{

	TABLE <- read.table(file, row.names = tr.col, header = header, skip = skip, blank.lines.skip = TRUE, comment.char = "#") 	

	if(!is.null(skip.col))
	{
		TABLE <- .RN_delete_col(TABLE, tr.col, skip.col)
	}

	if(!is.null(col.names))
	{
		colnames(TABLE) <- col.names
	}

	return (RN_calc(TABLE, design))
}
