
nullModelFastSKATTestPrep <- function(nullmod, threshold = 1e-10, idx.exclude = NULL){
    require(Matrix)
    
    if (is.null(idx.exclude)){
        Y <- nullmod$workingY
        W <- nullmod$model.matrix
        C <- nullmod$cholSigmaInv
    } else{
        Y <- nullmod$workingY[-idx.exclude]
        W <- nullmod$model.matrix[-idx.exclude,]
        C <- subsetCholSigmaInv(nullmod$cholSigmaInv, idx.exclude)
    }
    
    sparseC <- Matrix(C)	
	SIGMA <- chol2inv(t(C))
	SIGMA[abs(SIGMA) < threshold] <- 0
	SIGMA <- Matrix(SIGMA)
	CholSigma <- t(chol(SIGMA)) 
	SIGMAinv <- tcrossprod(sparseC)
	
	qr <- qr(as.matrix(solve(CholSigma,W)))


    out <- list(CholSigma = CholSigma, SIGMAinv = SIGMAinv, qr = qr, resids = nullmod$resid.marginal)
 #   class(out) <- "GENESIS.nullModelPrep"
    return(out)
}



.testVariantSetFastSKAT <- function(nullprep, G, weights){
    
    centeredG <- colMeans(G)
    ww <- weights(centeredG/2)
    sparseG <- Matrix(G, sparse = TRUE) %*% Diagonal(x = ww)
    
    rval <- list(mult = function(X){
    	base::qr.resid(nullprep$qr, as.matrix(solve(nullprep$CholSigma, (spG %*% X))))	
    }, tmult = function(X){
    	crossprod(spG, solve(t(nullprep$CholSigma), base::qr.resids(nullprep$qr,X)))
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
    
### or instead directly: ???
    stdres <- nullprep$SIGMAinv %*% nullprep$resids
	s = crossprod(sparseG, stdres)
	pval = pQF(sum(s^2), ...)  ## add arguments here??
}

