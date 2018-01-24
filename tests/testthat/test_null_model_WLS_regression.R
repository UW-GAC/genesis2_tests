context("check null model WLS regression")

.compareGENESIS <- function(gen1, gen2) {
    expect_equivalent(gen2$fixef, gen1$fixef, tolerance=1e-6)
    expect_equivalent(gen2$betaCov, gen1$betaCov, tolerance=1e-6)
    expect_equivalent(gen2$resid.marginal, gen1$resid.marginal, tolerance=1e-6)
    expect_equivalent(gen2$logLik, gen1$logLik, tolerance=1e-6)
    expect_equivalent(gen2$logLikR, gen1$logLikR, tolerance=1e-6)

## currently GENESIS has a mistake, in AUC calculation it uses the number of 
## matrices and groups used, but not the actual number of non-zero variance components. 
## so this "2" fixes for it. 
    expect_equivalent(gen2$AIC, gen1$AIC - 2, tolerance=1e-6)
    expect_equivalent(gen2$workingY, gen1$workingY, tolerance=1e-6)
    expect_equivalent(gen2$model.matrix, gen1$model.matrix, tolerance=1e-6)
    expect_equivalent(gen2$varComp, gen1$varComp[2:3], tolerance=1e-6)
    expect_equivalent(gen2$varCompCov, gen1$varCompCov[2:3, 2:3], tolerance=1e-4)
    expect_equivalent(gen2$family$family, gen1$family$family)
    expect_true(all(abs(gen2$cholSigmaInv - gen1$cholSigmaInv) < 1e-6))
    expect_equivalent(gen2$RSS, gen1$RSS, tolerance=1e-6)
}


test_that("WLS", {
### Checks for the linear regression case:
##### successful! 
n <- 100
X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))
y <- X %*% c(1, 0.5, 1) + rnorm(n, sd = c(rep(4, n/2), rep(2, n/2)))

group.idx <- list(G1 = c(1:(n/2)), G2 = c((n/2 + 1):n))
nullmod <- fitNullModel(y, X, group.idx = group.idx, verbose=FALSE)


## compare to GENESIS. I need to create a jusk correlation matrix that will zero out as having variance component zero...
scanData <- data.frame(scanID = paste0("p", 1:n), y = y, X1 = X[,1], X2 = X[,2], X3 = X[,3], group = c(rep("G1", n/2), rep("G2", n/2)))

varCompJunk <- TRUE
while(varCompJunk){
	cor.mat <- diag(rnorm(n, sd = 0.01))
	dimnames(cor.mat) <- list(scanData$scanID, scanData$scanID)
	lmm.genesis <- GENESIS::fitNullMM(scanData, "y", covars = c("X1", "X2", "X3"), covMatList = cor.mat,group.var = "group", verbose=FALSE)
	if (lmm.genesis$varComp[1] == 0 ) varCompJunk <- FALSE
}

expect_equal(nullmod$family$family, "gaussian")
expect_false(nullmod$family$mixedmodel)
expect_true(nullmod$hetResid)
expect_true(nullmod$converged)
expect_null(nullmod$zeroFLAG) 
expect_equivalent(nullmod$workingY, y)
expect_equivalent(nullmod$outcome, y)
expect_equivalent(nullmod$model.matrix, X)

## checks - GENESIS:
.compareGENESIS(lmm.genesis, nullmod)

})
