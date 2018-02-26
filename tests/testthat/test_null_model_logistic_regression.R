context("check null model logistic regression")

test_that("logistic", {
### Checks for the logistic regression case:
##### successful! 
n <- 100
X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))

expit <- function(x){exp(x)/(1+exp(x))}
p <- expit(X %*% c(-1, 0.5, 1))
D <- rbinom(n, size = 1, prob = p)

nullmod <- fitNullMod(D, X, family = "binomial", verbose=FALSE)

# compare to logistic regression
glm.mod <- glm(D ~ -1 + X, family = "binomial")

## checks - logistic regression:
expect_equal(nullmod$family$family, "binomial")
expect_false(nullmod$hetResid)
expect_equivalent(nullmod$fitted.values, fitted(glm.mod))
expect_false(nullmod$family$mixedmodel)
expect_equivalent(nullmod$resid.marginal, resid(glm.mod, type = "response"))
expect_true(all(nullmod$fixef == summary(glm.mod)$coef))
expect_equivalent(nullmod$varComp, fitted(glm.mod)*(1-fitted(glm.mod)))
expect_null(nullmod$varCompCov)
expect_equivalent(nullmod$betaCov, vcov(glm.mod))
expect_equivalent(nullmod$fitted.values, fitted(glm.mod))
expect_equal(nullmod$logLik, as.numeric(logLik(glm.mod)))
expect_equal(nullmod$AIC, AIC(glm.mod))
expect_equivalent(nullmod$workingY, D)
expect_equivalent(nullmod$outcome, D)
expect_equivalent(nullmod$model.matrix, X)
expect_equivalent(nullmod$cholSigmaInv, 1/sqrt(fitted(glm.mod)*(1-fitted(glm.mod))))
expect_equal(nullmod$converged, glm.mod$converged)
expect_null(nullmod$zeroFLAG)
expect_equal(nullmod$RSS, 1)


## compare to GENESIS:
scanData <- data.frame(scanID = paste0("p", 1:n), D = D, X1 = X[,1], X2 = X[,2], X3 = X[,3])
glm.genesis <- GENESIS::fitNullReg(scanData, "D", covars = c("X1", "X2", "X3"), family = "binomial", verbose=FALSE)


## checks - GENESIS:
expect_equivalent(nullmod$fixef, glm.genesis$fixef)
expect_equivalent(nullmod$betaCov, glm.genesis$betaCov)
expect_equivalent(nullmod$resid.marginal, glm.genesis$resid.response)
expect_equivalent(nullmod$logLik, glm.genesis$logLik)
expect_equivalent(nullmod$AIC, glm.genesis$AIC)
expect_equivalent(nullmod$workingY, glm.genesis$workingY)
expect_equivalent(nullmod$model.matrix, glm.genesis$model.matrix)
expect_equal(nullmod$family$family, glm.genesis$family$family)
expect_equivalent(nullmod$varComp, glm.genesis$sigma^2)

})
