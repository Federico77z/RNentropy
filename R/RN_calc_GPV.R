RN_calc_GPV <-
function (X, bind=TRUE) 
{
    	rnums <- sapply(X, is.numeric)
    	GL_LPV <- apply(X[rnums], 1, '.RN_calc_GPV_row')
    	if(bind)
    	{
    		GPV <- cbind(X, GL_LPV)
    		return(GPV)
	}
	else
	{
		return(GL_LPV)
	}
}
