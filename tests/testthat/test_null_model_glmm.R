context("check null model glmm")

.compareGENESIS <- function(gen1, gen2) {
    expect_equivalent(gen2$fixef, gen1$fixef, tolerance=1e-6)
    expect_equivalent(gen2$betaCov, gen1$betaCov, tolerance=1e-6)
    expect_equivalent(gen2$resid.marginal, gen1$resid.marginal, tolerance=1e-6)
    expect_equivalent(gen2$logLik, gen1$logLik, tolerance=1e-6)
    expect_equivalent(gen2$logLikR, gen1$logLikR, tolerance=1e-6)

    expect_equivalent(gen2$AIC, gen1$AIC, tolerance=1e-6)
    expect_equivalent(gen2$workingY, gen1$workingY, tolerance=1e-6)
    expect_equivalent(gen2$model.matrix, gen1$model.matrix, tolerance=1e-6)
    expect_equivalent(gen2$varComp, gen1$varComp, tolerance=1e-6)
    expect_equivalent(gen2$varCompCov, gen1$varCompCov, tolerance=1e-4)  ## not smaller then 1e-6 because values were not updated after convergence. However, this is not important, because the variance component values are converged at the desired level!
    expect_equivalent(gen2$family$family, gen1$family$family)
    expect_equivalent(gen2$zeroFLAG, gen1$zeroFLAG)
    expect_true(all(abs(gen2$cholSigmaInv - gen1$cholSigmaInv) < 1e-6))
}


test_that("glmm", {
### Checks for the logistic (also poisson, same algorithm) regression case:
n <- 100
scanID <- paste0("p", 1:n)
X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))

## generate random effects:
sqrt.cor.mat <- matrix(rnorm(n*n, sd = 0.05),n,n)
cor.mat <- crossprod(sqrt.cor.mat)
dimnames(cor.mat) <- list(scanID, scanID)


## compare to GENESIS. Make sure the variance component isn't zeroed out.
## if it does, the function has an error. 

varCompZero <- TRUE
while(varCompZero){	
	random.iid <- rnorm(n)
	random <- crossprod(sqrt.cor.mat*0.05, random.iid)
	expit <- function(x){exp(x)/(1+exp(x))} 
	p <- expit(X %*% c(-1, 0.5, 1) + random) 
	D <- rbinom(n, size = 1, prob = p)
	
	scanData <- data.frame(scanID = scanID, D = D, X1 = X[,1], X2 = X[,2], X3 = X[,3])
	
	glmm.genesis <- tryCatch({
		GENESIS::fitNullMM(scanData, "D", covars = c("X1", "X2", "X3"), covMatList = cor.mat, family = "binomial", verbose=FALSE)
				}, 
			warning = function(w){return(list(message = "warning"))},
			error = function(e){return(list(message = "error"))}
			)
	if (!is.null(glmm.genesis$message)) next
	if (glmm.genesis$varComp[1] != 0 ) varCompZero <- FALSE
}

nullmod <- fitNullMod(y = D, X = X, covMatList = cor.mat, verbose=FALSE, family = "binomial")

expect_equal(nullmod$family$family, "binomial")
expect_true(nullmod$family$mixedmodel)
expect_false(nullmod$hetResid)
expect_true(nullmod$converged)
expect_false(nullmod$zeroFLAG)
expect_equivalent(nullmod$outcome, D)
expect_equivalent(nullmod$model.matrix, X)
expect_equal(nullmod$RSS, 1)

## checks - GENESIS:
.compareGENESIS(glmm.genesis, nullmod)

})
