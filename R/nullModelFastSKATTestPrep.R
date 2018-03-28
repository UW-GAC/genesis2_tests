
nullModelFastSKATTestPrep <- function(nullmod, threshold = 1e-10){
    
    Y <- nullmod$workingY
    W <- nullmod$model.matrix
    C <- nullmod$cholSigmaInv

	## C can also be a scalar (simple linear regression) - adapt code  
	if (is.null(dim(C)) & length(C) == 1){
		CholSigma <- 1/C
		SIGMAinv <- C^2
		qr <- qr(W/CholSigma)
	} else{
		
  	  SIGMA <- chol2inv(t(C))
  	  SIGMA[abs(SIGMA) < threshold] <- 0
   	 SIGMA <- Matrix(SIGMA)
   	 CholSigma <- t(chol(SIGMA)) 
   	 SIGMAinv <- tcrossprod(C)

   	 qr <- qr(as.matrix(solve(CholSigma,W)))
	}

    out <- list(CholSigma = CholSigma, SIGMAinv = SIGMAinv, qr = qr, resids = nullmod$resid.marginal)
    return(out)
}



.testVariantSetFastSKAT <- function(nullprep, G, weights,  rho = 0, method=c("ssvd","lanczos","satterthwaite"), neig=100, tr2.sample.size=500, q=NULL,  convolution.method=c("saddlepoint","integration"), remainder.underflow=c("warn","missing","error")){
    
    if (rho == 0){
		sparseG <- Matrix(G, sparse = TRUE) %*% Diagonal(x = weights) 
	} else{
		rho.mat <- matrix(rho, nrow = length(weights), ncol = length(weights))
		diag(rho.mat) <- 1
		cholRhoMat <- t(chol(rho.mat, pivot=TRUE))
		sparseG <- Matrix(G, sparse = TRUE) %*% Diagonal(x = weights)  %*% cholRhoMat
	}
    
	if (is.null(dim(nullprep$CholSigma)) & length(nullprep$CholSigma) == 1){
		  rval <- list(mult = function(X){
    		base::qr.resid(nullprep$qr, as.matrix((sparseG %*% X)/nullprep$CholSigma))	
	    }, tmult = function(X){
	    	crossprod(sparseG, base::qr.resid(nullprep$qr,X)/nullprep$CholSigma)
	    }, 
	    trace = NULL,
	    ncol = ncol(G),
	    nrow = nrow(G),
	    Q = function(){
	    	stdres <- nullprep$SIGMAinv * nullprep$resids
	    	s = crossprod(sparseG, stdres)
	    	sum(s^2)
	    })
	    
	    class(rval) <- c( "matrixfree") ## or famSKAT_genesis?
	
		} else{
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

			class(rval) <- c("famSKAT_genesis", "famSKAT", "matrixfree") 
			
			}
    

    
    pval <- bigQF::pQF(rval$Q(), rval,  method = method, neig = neig, tr2.sample.size = tr2.sample.size, q = q, convolution.method = convolution.method, remainder.underflow = remainder.underflow) ### what object should be the second entry of pQF??
   	names(pval) <- "pval"
   	
	return(pval)
}




