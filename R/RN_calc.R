RN_calc <-
function(X, design = NULL)
{
  if(is.null(design))
  {
    design <- .RN_default_design(sum(sapply(X, is.numeric)))
  }
  
	Results <- list(expr = X, design = design)

	GPV <- RN_calc_GPV(X, bind = FALSE)
	LPV <- RN_calc_LPV(X, design = design, bind = FALSE)	

	TABLE = cbind(X,'---',GPV,'---',LPV)	

	Results$gpv <- GPV
	Results$lpv <- LPV
	Results$c_like <- TABLE
	Results$res <- cbind(X, GPV, LPV)
	
	return(Results)
}
