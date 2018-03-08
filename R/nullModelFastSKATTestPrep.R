
nullModelFastSKATTestPrep <- function(nullmod, threshold = 1e-10){
    
    Y <- nullmod$workingY
    W <- nullmod$model.matrix
    C <- nullmod$cholSigmaInv
    
    #sparseC <- Matrix(C) # should be already a Matrix
    SIGMA <- chol2inv(t(C))
    SIGMA[abs(SIGMA) < threshold] <- 0
    SIGMA <- Matrix(SIGMA)
    CholSigma <- t(chol(SIGMA)) 
    SIGMAinv <- tcrossprod(C)
    
    qr <- qr(as.matrix(solve(CholSigma,W)))


    out <- list(CholSigma = CholSigma, SIGMAinv = SIGMAinv, qr = qr, resids = nullmod$resid.marginal)
    return(out)
}



.testVariantSetFastSKAT <- function(nullprep, G, weights){
    
    sparseG <- Matrix(G, sparse = TRUE) %*% Diagonal(x = weights)
    
    rval <- list(mult = function(X){
    	base::qr.resid(nullprep$qr, as.matrix(solve(nullprep$CholSigma, (sparseG %*% X))))	
    }, tmult = function(X){
    	crossprod(sparseG, solve(t(nullprep$CholSigma), base::qr.resid(nullprep$qr,X)))
    }, 
    trace = NULL,
    ncol = ncol(G),
    nrow = nrow(G),
    Q = function(){
    	stdres <- nullprep$SIGMAinv %*% nullprep$resids
    	s = crossprod(sparseG, stdres)
    	sum(s^2)
    })
    
    class(rval) <- c("famSKAT_lmekin", "famSKAT", "matrixfree") ## or famSKAT_genesis?
    return(rval)
    
    pval <- bigQF::pQF(rval$Q(), rval, neig = ) ### what object should be the second entry of pQF??
    
    
# ### or instead directly: ???
    # stdres <- nullprep$SIGMAinv %*% nullprep$resids
	# s = crossprod(sparseG, stdres)
	# pval = pQF(sum(s^2), ...)  ## add arguments here??
}



