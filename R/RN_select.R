RN_select <-
function(Results, gpv_t = 0.01, lpv_t = 0.01, method = "BH")
{
	lpv_t <- -log10(lpv_t)
	gpv_t <- -log10(gpv_t)

	Results$gpv_bh <- -log10(p.adjust(10 ^ -Results$gpv, method=method))
	
	true_rows <- (Results$gpv_bh >= gpv_t)

	design_b <- t(Results$design > 0)	
	
	Results$lpv_sel <- data.frame(row.names = rownames(Results$lpv)[true_rows])
	
	boh <- list()
	
	for(d in seq_along(design_b[,1]))
	{
	  Results$lpv_sel <- cbind(Results$lpv_sel, 
	                           apply(as.matrix(Results$lpv[true_rows,design_b[d,, drop = F]]), 
	                                 1, ".RN_select_lpv_row", lpv_t))
	}
	
	colnames(Results$lpv_sel) <- colnames(Results$design)

	lbl <- Results$res[,!sapply(Results$res, is.numeric), drop = FALSE]

	Results$selected <- cbind(lbl[true_rows,], Results$gpv[true_rows], Results$gpv_bh[true_rows], Results$lpv_sel)

	colnames(Results$selected) <- c(names(which(!sapply(Results$res, is.numeric))), "GL_LPV", "Corr. GL_LPV", colnames(Results$lpv_sel))

	Results$selected <- Results$selected[order(Results$selected[,3], decreasing=TRUE),]

	Results$lpv_sel <- NULL
	
	return (Results)
}
