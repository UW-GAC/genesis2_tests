context("check update of null model after ranknormalizing and re-scaling residuals")
require(GENESIS)
require(GWASTools)

test_that("updateOutcome", {
### Checks that updating the outcome and re-fitting the null model works okay.
n <- 100
X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))
y <- X %*% c(1, 0.5, 1) + rnorm(n, sd = c(rep(4, n/2), rep(2, n/2)))

scanID <- paste0("p", 1:n)
group.idx <- list(G1 = c(1:(n/2)), G2 = c((n/2 + 1):n))

cor.mat <- matrix(rnorm(n*n, sd = 0.05),n,n)
cor.mat <- crossprod(cor.mat)
dimnames(cor.mat) <- list(scanID, scanID)
covMatList <- list(A = cor.mat)

nullmod <- fitNullModel(y, X, group.idx = group.idx, covMatList, verbose=FALSE)

group.ind <- 1
expect_equal(.averageGroupVar(nullmod$varComp, covMatList = covMatList, group.idx = group.idx[group.ind]),
             nullmod$varComp[1]*mean(diag(cor.mat)[group.idx[[group.ind]]]) + nullmod$varComp[group.ind + 1])

group.ind <- 2
expect_equal(.averageGroupVar(nullmod$varComp, covMatList = covMatList, group.idx = group.idx[group.ind]),
	     nullmod$varComp[1]*mean(diag(cor.mat)[group.idx[[group.ind]]]) + nullmod$varComp[group.ind + 1])

expect_equal(.averageGroupVar(nullmod$varComp, covMatList = NULL, group.idx = group.idx[group.ind]),
             nullmod$varComp[group.ind + 1])
						
throws_error(.averageGroupVar(nullmod$varComp, covMatList = NULL, group.idx = group.idx))
throws_error(.averageGroupVar(nullmod$varComp, covMatList = NULL, group.idx = NULL))
throws_error(.averageGroupVar(nullmod$varComp, covMatList = NULL, group.idx = group.idx[[1]]))

nullmod2 <- updateNullModOutcome(nullmod, covMatList = covMatList, group.idx = group.idx,  rankNorm.option = c("by.group"), rescale = c("None"), verbose=FALSE)


expect_true(abs(.averageGroupVar(nullmod2$varComp, covMatList = covMatList, group.idx = group.idx[1]) - 1) < 0.1 )
expect_true(abs(.averageGroupVar(nullmod2$varComp, covMatList = covMatList, group.idx = group.idx[2]) -1) < 0.1 )

# some comparisons when the residuals are rank-normalized and re-scaled within gorups
average.group1.model.var.nullmod <- .averageGroupVar(nullmod$varComp, covMatList = covMatList, group.idx = group.idx[1])
average.group2.model.var.nullmod <- .averageGroupVar(nullmod$varComp, covMatList = covMatList, group.idx = group.idx[2])

nullmod3 <- updateNullModOutcome(nullmod, covMatList = covMatList, group.idx = group.idx,  rankNorm.option = c("by.group"), rescale = c("residSD"), verbose=FALSE)
nullmod4 <- updateNullModOutcome(nullmod, covMatList = covMatList, group.idx = group.idx,  rankNorm.option = c("by.group"), rescale = c("model"), verbose=FALSE)

average.group1.model.var.nullmod3 <- .averageGroupVar(nullmod3$varComp, covMatList = covMatList, group.idx = group.idx[1])
average.group2.model.var.nullmod3 <- .averageGroupVar(nullmod3$varComp, covMatList = covMatList, group.idx = group.idx[2])
average.group1.model.var.nullmod4 <- .averageGroupVar(nullmod4$varComp, covMatList = covMatList, group.idx = group.idx[1])
average.group2.model.var.nullmod4 <- .averageGroupVar(nullmod4$varComp, covMatList = covMatList, group.idx = group.idx[2])

# compare to model 3:
expect_true(average.group1.model.var.nullmod3/average.group1.model.var.nullmod < 1.1 & average.group1.model.var.nullmod3/average.group1.model.var.nullmod > 0.9 )
expect_true(average.group2.model.var.nullmod3/average.group2.model.var.nullmod < 1.1 & average.group2.model.var.nullmod3/average.group2.model.var.nullmod > 0.9 )

# compare to model 4:
expect_true(average.group1.model.var.nullmod4/average.group1.model.var.nullmod < 1.1 & average.group1.model.var.nullmod4/average.group1.model.var.nullmod > 0.9 )
expect_true(average.group2.model.var.nullmod4/average.group2.model.var.nullmod < 1.1 & average.group2.model.var.nullmod4/average.group2.model.var.nullmod > 0.9 )



})
