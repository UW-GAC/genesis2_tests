
## calculate BIC given a parameter vector beta, assuming that all variance components are fixed. 
## the BIC is simplifed because we don't need to compare (and therefore to calculate) fixed values between models.
## we will need to decide whether this function takes nullmod, or arguments of nullmod.
## perhaps call the function .calcBICcompareModels or something that reflects that his is not really the BIC,
## but rather only an additive value that changes between models under comparison?

.calcBICfixedVCs <- function(nullmod, G, beta){
	n <- length(nullmod$workingY)
	k <- length(beta)  + nrow(nullmod$fixef) 
	
	Sigma.inv <- crossprod(nullmod$cholSigmaInv)

    ## calculate the mean of the outcomes
    fits <- tcrossprod(X, t(nullmod$fixef[,1])) + tcrossprod(G, t(beta))
    
    # obtain marginal residuals
    residM <- as.vector(Y - fits)
    Sigma.inv_R <- crossprod(Sigma.inv, residM)
    
    # calculate the sum of squares part of the likelihood
    Rt_Sigma.inv_R <- crossprod(residM, Sigma.inv_R)
    
    ## I think we only need this part, given the assumption on fixed variance components: 
    BIC <- log(n)*k + Rt_Sigma.inv_R
    
    return(BIC)
    
    #RSS <- as.numeric(Rt_Sigma.inv_R/(n - k)) 
    #logLik <- as.numeric(-0.5 * n * log(2 * pi * RSS) - 0.5 * Rt_Sigma.inv_R/RSS)
	
	
	
}