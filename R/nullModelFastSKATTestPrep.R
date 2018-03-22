
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



.testVariantSetFastSKAT <- function(nullprep, G, weights, method=c("ssvd","lanczos","satterthwaite"), neig=100, tr2.sample.size=500, q=NULL,  convolution.method=c("saddlepoint","integration"), remainder.underflow=c("warn","missing","error")){
    
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
    
    class(rval) <- c("famSKAT_genesis", "famSKAT", "matrixfree") ## or famSKAT_genesis?
    
    pval <- bigQF::pQF(rval$Q(), rval,  method = method, neig = neig, tr2.sample.size = tr2.sample.size, q = q, convolution.method = convolution.method, remainder.underflow = remainder.underflow) ### what object should be the second entry of pQF??
   	names(pval) <- "pval"
   	
	return(pval)
}




