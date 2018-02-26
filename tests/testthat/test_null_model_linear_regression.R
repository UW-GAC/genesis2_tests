context("check null model linear regression")

test_that("linear", {
### Checks for the linear regression case:
##### successful! 
n <- 100
X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))
y <- X %*% c(1, 0.5, 1) + rnorm(n, sd = 2)

nullmod <- fitNullMod(y, X, verbose=FALSE)

# compare to linear regression
lm.mod <- lm(y ~ -1 + X)

## checks - linear regression:
expect_equal(nullmod$fitted.values, fitted(lm.mod))
expect_equal(nullmod$family$family, "gaussian")
expect_false(nullmod$family$mixedmodel)
expect_false(nullmod$hetResid)
expect_equivalent(nullmod$resid.marginal, lm.mod$resid)
expect_true(all(nullmod$fixef == summary(lm.mod)$coef))
expect_equal(nullmod$varComp, summary(lm.mod)$sigma^2)
expect_null(nullmod$varCompCov)
expect_equivalent(nullmod$betaCov, vcov(lm.mod))
expect_equivalent(nullmod$fitted.values, fitted(lm.mod))
expect_equal(nullmod$logLik, as.numeric(logLik(lm.mod)))
expect_equal(nullmod$AIC, AIC(lm.mod))
expect_equivalent(nullmod$workingY, y)
expect_equivalent(nullmod$outcome, y)
expect_equivalent(nullmod$model.matrix, X)
expect_equal(nullmod$cholSigmaInv, 1/summary(lm.mod)$sigma)
expect_true(nullmod$converged)
expect_null(nullmod$zeroFLAG)
expect_equal(nullmod$RSS, sum(lm.mod$resid^2)/(summary(lm.mod)$sigma^2*(n - ncol(X))))

## compare to GENESIS:
scanData <- data.frame(scanID = paste0("p", 1:n), y = y, X1 = X[,1], X2 = X[,2], X3 = X[,3])
lm.genesis <- GENESIS::fitNullReg(scanData, "y", covars = c("X1", "X2", "X3"), verbose=FALSE)


## checks - GENESIS:
expect_equivalent(nullmod$fixef, lm.genesis$fixef)
expect_equivalent(nullmod$betaCov, lm.genesis$betaCov)
expect_equivalent(nullmod$resid.marginal, lm.genesis$resid.response)
expect_equivalent(nullmod$logLik, lm.genesis$logLik)
expect_equivalent(nullmod$AIC, lm.genesis$AIC)
expect_equivalent(nullmod$workingY, lm.genesis$workingY)
expect_equivalent(nullmod$model.matrix, lm.genesis$model.matrix)
expect_equivalent(nullmod$varComp, lm.genesis$sigma^2)
expect_equal(nullmod$family$family, lm.genesis$family$family)
})

