context("check null model lmm")

.compareGENESIS <- function(gen1, gen2) {
    expect_equivalent(gen2$fixef, gen1$fixef)
    expect_equivalent(gen2$betaCov, gen1$betaCov)
    expect_equivalent(gen2$resid.marginal, gen1$resid.marginal)
    expect_equivalent(gen2$logLik, gen1$logLik)
    expect_equivalent(gen2$logLikR, gen1$logLikR)

    expect_equivalent(gen2$AIC, gen1$AIC)
    expect_equivalent(gen2$workingY, gen1$workingY)
    expect_equivalent(gen2$model.matrix, gen1$model.matrix)
    expect_equivalent(gen2$varComp, gen1$varComp)
    expect_equivalent(gen2$varCompCov, gen1$varCompCov)
    expect_equivalent(gen2$family$family, gen1$family$family)
    expect_equivalent(gen2$zeroFLAG, gen1$zeroFLAG)
    expect_true(all(abs(gen2$cholSigmaInv - gen1$cholSigmaInv) < 1e-9))
    expect_equivalent(gen2$RSS, gen1$RSS)
}

test_that("lmm", {
### Checks for the linear regression case:
##### successful! 
n <- 100
X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))
y <- X %*% c(1, 0.5, 1) + rnorm(n, sd = c(rep(4, n/2), rep(2, n/2)))

group.idx <- list(G1 = c(1:(n/2)), G2 = c((n/2 + 1):n))

## compare to GENESIS. I need to create a jusk correlation matrix that will zero out as having variance component zero...
scanData <- data.frame(scanID = paste0("p", 1:n), y = y, X1 = X[,1], X2 = X[,2], X3 = X[,3], group = c(rep("G1", n/2), rep("G2", n/2)))

## varCompJunk <- FALSE
## while(!varCompJunk){
	cor.mat <- matrix(rnorm(n*n, sd = 0.05),n,n)
	cor.mat <- crossprod(cor.mat)
	dimnames(cor.mat) <- list(scanData$scanID, scanData$scanID)
	lmm.genesis <- GENESIS::fitNullMM(scanData, "y", covars = c("X1", "X2", "X3"), covMatList = cor.mat,group.var = "group", verbose=FALSE)
## 	if (lmm.genesis$varComp[1] != 0 ) varCompJunk <- TRUE
## }

nullmod <- fitNullMod(y, X, group.idx = group.idx, cor.mat, verbose=FALSE)


expect_equal(nullmod$family$family, "gaussian")
expect_true(nullmod$family$mixedmodel)
expect_true(nullmod$hetResid)
expect_true(nullmod$converged)
expect_equivalent(nullmod$workingY, y)
expect_equivalent(nullmod$outcome, y)
expect_equivalent(nullmod$model.matrix, X)


## checks - GENESIS:
.compareGENESIS(lmm.genesis, nullmod)


### test updating a conditional model: 
G = matrix(rnorm(100, 100,1))
nullmod2 <- updateNullModCond(nullmod, G, covMatList = list(cor.mat), AIREML.tol = 1e-7, verbose=FALSE)
nullmod3 <- fitNullMod(y, cbind(X, G), group.idx = group.idx, cor.mat, AIREML.tol = 1e-7, verbose=FALSE)

expect_equivalent(nullmod2$varComp, nullmod3$varComp, tolerance=1e-5)
expect_equivalent(nullmod3$fixef, nullmod2$fixef, tolerance=1e-5)
expect_equivalent(nullmod2$cholSigmaInv, nullmod3$cholSigmaInv, tolerance=1e-5)
expect_equivalent(nullmod2$varCompCov, nullmod3$varCompCov, tolerance=1e-5)




## test without group
nullmod <- fitNullMod(y, X, cor.mat, verbose=FALSE)

lmm.genesis <- GENESIS::fitNullMM(scanData, "y", covars = c("X1", "X2", "X3"), covMatList = cor.mat, verbose=FALSE)

expect_false(nullmod$hetResid)
expect_true(nullmod$converged)
expect_true(all(nullmod$workingY == y))
expect_true(all(nullmod$outcome == y))
expect_equivalent(nullmod$workingY, y)
expect_equivalent(nullmod$outcome, y)
expect_equivalent(nullmod$model.matrix, X)


## checks - GENESIS:
.compareGENESIS(lmm.genesis, nullmod)

})
