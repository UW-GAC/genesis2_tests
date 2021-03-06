context("check variant set association tests")
library(bigQF)

.prepNullmod <- function(n, MM=FALSE, binary=FALSE, het.var = TRUE) {
    X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))
    if (!binary) {
	if (het.var) {
		y <- X %*% c(1, 0.5, 1) + rnorm(n, sd = c(rep(4, n/2), rep(2, n/2)))
		} else{
			y <- X %*% c(1, 0.5, 1) + rnorm(n, sd = rep(4, n))
		}
        family <- "gaussian"
    }  else {
        y <- rbinom(n, size = 1, prob = 0.4)
        family <- "binomial"
    }
	
	group.idx <- list(G1 = c(1:(n/2)), G2 = c((n/2 + 1):n))

        if (MM) {
            cor.mat <- matrix(rnorm(n*n, sd = 0.05),n,n)
            cor.mat <- crossprod(cor.mat)
            covMatList <- list(A = cor.mat)
        } else {
            covMatList <- NULL
        }
	if (het.var){
		nullmod <- fitNullMod(y, X, covMatList = covMatList, group.idx = group.idx, family=family, verbose=FALSE)
	} else{
		nullmod <- fitNullMod(y, X, covMatList = covMatList, family=family, verbose=FALSE)
	}
	return(nullmod)
	
	#
}


test_that("SKAT with rho=1 matches burden", {
	n <- 100

	### create a matrix of genetic variants to test.
	geno <- matrix(rbinom(200*n, size = 2, prob = 0.2), nrow = n, ncol = 200)
        weights <- rep(1, ncol(geno))

        ## mixed model
        nullmod <- .prepNullmod(n, MM=TRUE)
        skat.1 <- .testVariantSetSKAT(nullmod, geno, weights, rho=0, pval.method="davies")
       
       nullprep.2 <- nullModelFastSKATTestPrep(nullmod)
       skat.2 <- .testVariantSetFastSKAT(nullprep.2, geno, weights) 
        
       expect_true(abs(skat.1$pval_0 - skat.2) < 0.01)
        
        ## basic - WLS
        nullmod <- .prepNullmod(n, MM=FALSE)
        skat.1 <- .testVariantSetSKAT(nullmod, geno, weights, rho=0, pval.method="davies")
       
       nullprep.2 <- nullModelFastSKATTestPrep(nullmod)
       skat.2 <- .testVariantSetFastSKAT(nullprep.2, geno, weights) 
        
       expect_true(abs(skat.1$pval_0 - skat.2) < 0.01)
        
        ## basic linear reg
        nullmod <- .prepNullmod(n, MM=FALSE, het.var = FALSE)
        skat.1 <- .testVariantSetSKAT(nullmod, geno, weights, rho=0, pval.method="davies")
       
       nullprep.2 <- nullModelFastSKATTestPrep(nullmod)
       skat.2 <- .testVariantSetFastSKAT(nullprep.2, geno, weights) 
        
       expect_true(abs(skat.1$pval_0 - skat.2) < 0.01)


        
        ## mixed model - binary
	 nullmod <- .prepNullmod(n, MM=TRUE, binary=TRUE)
        skat.1 <- .testVariantSetSKAT(nullmod, geno, weights, rho=0, pval.method="davies")
       
       nullprep.2 <- nullModelFastSKATTestPrep(nullmod)
       skat.2 <- .testVariantSetFastSKAT(nullprep.2, geno, weights) 
        
       expect_true(abs(skat.1$pval_0 - skat.2) < 0.01)

        
        ## basic - binary
	nullmod <- .prepNullmod(n, MM=FALSE, binary=TRUE)
        skat.1 <- .testVariantSetSKAT(nullmod, geno, weights, rho=0, pval.method="davies")
       
       nullprep.2 <- nullModelFastSKATTestPrep(nullmod)
       skat.2 <- .testVariantSetFastSKAT(nullprep.2, geno, weights) 
        
       expect_true(abs(skat.1$pval_0 - skat.2) < 0.01)
       
       
       ### now incorporate rho > 0:
      nullmod <- .prepNullmod(n, MM=TRUE)
        skat.1 <- .testVariantSetSKAT(nullmod, geno, weights, rho=0.3, pval.method="davies")  
        nullprep.2 <- nullModelFastSKATTestPrep(nullmod)
       skat.2 <- .testVariantSetFastSKAT(nullprep.2, geno, weights, rho = 0.3) 

		 expect_true(abs(skat.1$pval_0 - skat.2) < 0.01)

})
