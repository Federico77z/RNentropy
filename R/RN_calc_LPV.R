RN_calc_LPV <-
function(X, design = NULL, bind = TRUE)
{
  if(is.null(design))
  {
    design <- .RN_default_design(sum(sapply(X, is.numeric)))
  }
  
	rnums <- sapply(X, is.numeric)

	.RN_design_check(X, design, rnums)
	RL <- .RN_get_replicate_list(design)
	
	l1 <- rep("LOC_LPV", length(X[rnums]))
	l2 <- colnames(X[rnums])
	PV <- apply(X[rnums],  1, '.RN_calc_LPV_row', RL = RL)
	PV <- t(PV)
	colnames(PV) <- paste(l1, l2, sep="_")

	if(bind)
	{
		return (cbind(X, PV))
	}
	else
		return (PV)
}
