context("compare to GENESIS")

test_that("fitNullMod matches fitNullReg - linear", {
    dat <- .testNullInputs()
    nullmod <- fitNullMod(dat$y, dat$X, verbose=FALSE)

    df <- .asDataFrame(dat)
    lm.genesis <- GENESIS::fitNullReg(df, outcome="y", covars=c("X1", "X2", "X3"), verbose=FALSE)

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


test_that("fitNullMod matches fitNullReg - binary", {
    dat <- .testNullInputs(binary=TRUE)
    nullmod <- fitNullMod(dat$y, dat$X, family="binomial", verbose=FALSE)

    df <- .asDataFrame(dat)
    glm.genesis <- GENESIS::fitNullReg(df, outcome="y", covars=c("X1", "X2", "X3"), family="binomial", verbose=FALSE)
    
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


test_that("fitNullMod matches fitNullMM - linear, no group", {
    dat <- .testNullInputs()
    nullmod <- fitNullMod(dat$y, dat$X, dat$cor.mat, verbose=FALSE)

    df <- .asDataFrame(dat)
    lmm.genesis <- GENESIS::fitNullMM(df, outcome="y", covars=c("X1", "X2", "X3"), covMatList=dat$cor.mat, verbose=FALSE)

    expect_equivalent(nullmod$fixef, lmm.genesis$fixef)
    expect_equivalent(nullmod$betaCov, lmm.genesis$betaCov)
    expect_equivalent(nullmod$resid.marginal, lmm.genesis$resid.marginal)
    expect_equivalent(nullmod$logLik, lmm.genesis$logLik)
    expect_equivalent(nullmod$logLikR, lmm.genesis$logLikR)

    expect_equivalent(nullmod$AIC, lmm.genesis$AIC)
    expect_equivalent(nullmod$workingY, lmm.genesis$workingY)
    expect_equivalent(nullmod$model.matrix, lmm.genesis$model.matrix)
    expect_equivalent(nullmod$varComp, lmm.genesis$varComp)
    expect_equivalent(nullmod$varCompCov, lmm.genesis$varCompCov)
    expect_equivalent(nullmod$family$family, lmm.genesis$family$family)
    expect_equivalent(nullmod$zeroFLAG, lmm.genesis$zeroFLAG)
    expect_true(all(abs(nullmod$cholSigmaInv - lmm.genesis$cholSigmaInv) < 1e-9))
    expect_equivalent(nullmod$RSS, lmm.genesis$RSS)
})


test_that("fitNullMod matches fitNullMM - linear, with group", {
    dat <- .testNullInputs()
    nullmod <- fitNullMod(dat$y, dat$X, dat$cor.mat, group.idx=dat$group.idx, verbose=FALSE)

    df <- .asDataFrame(dat)
    lmm.genesis <- GENESIS::fitNullMM(df, outcome="y", covars=c("X1", "X2", "X3"), covMatList=dat$cor.mat, group.var="group", verbose=FALSE)

    expect_equivalent(nullmod$fixef, lmm.genesis$fixef)
    expect_equivalent(nullmod$betaCov, lmm.genesis$betaCov)
    expect_equivalent(nullmod$resid.marginal, lmm.genesis$resid.marginal)
    expect_equivalent(nullmod$logLik, lmm.genesis$logLik)
    expect_equivalent(nullmod$logLikR, lmm.genesis$logLikR)

    expect_equivalent(nullmod$AIC, lmm.genesis$AIC)
    expect_equivalent(nullmod$workingY, lmm.genesis$workingY)
    expect_equivalent(nullmod$model.matrix, lmm.genesis$model.matrix)
    expect_equivalent(nullmod$varComp, lmm.genesis$varComp)
    expect_equivalent(nullmod$varCompCov, lmm.genesis$varCompCov)
    expect_equivalent(nullmod$family$family, lmm.genesis$family$family)
    expect_equivalent(nullmod$zeroFLAG, lmm.genesis$zeroFLAG)
    expect_true(all(abs(nullmod$cholSigmaInv - lmm.genesis$cholSigmaInv) < 1e-9))
    expect_equivalent(nullmod$RSS, lmm.genesis$RSS)
})


test_that("fitNullMod matches fitNullMM - binary", {
    varCompZero <- TRUE
    while(varCompZero){
        dat <- .testNullInputs(binary=TRUE)
        df <- .asDataFrame(dat)
        
        glmm.genesis <- tryCatch({
            GENESIS::fitNullMM(df, outcome="y", covars = c("X1", "X2", "X3"), covMatList=dat$cor.mat, family="binomial", verbose=FALSE)
        }, 
        warning = function(w){return(list(message = "warning"))},
        error = function(e){return(list(message = "error"))}
        )
        if (!is.null(glmm.genesis$message)) next
        if (glmm.genesis$varComp[1] != 0 ) varCompZero <- FALSE
    }

    nullmod <- fitNullMod(dat$y, dat$X, covMatList=dat$cor.mat, family="binomial", verbose=FALSE)
    
    expect_equivalent(nullmod$fixef, glmm.genesis$fixef, tolerance=1e-6)
    expect_equivalent(nullmod$betaCov, glmm.genesis$betaCov, tolerance=1e-6)
    expect_equivalent(nullmod$resid.marginal, glmm.genesis$resid.marginal, tolerance=1e-6)
    expect_equivalent(nullmod$logLik, glmm.genesis$logLik, tolerance=1e-6)
    expect_equivalent(nullmod$logLikR, glmm.genesis$logLikR, tolerance=1e-6)

    expect_equivalent(nullmod$AIC, glmm.genesis$AIC, tolerance=1e-6)
    expect_equivalent(nullmod$workingY, glmm.genesis$workingY, tolerance=1e-6)
    expect_equivalent(nullmod$model.matrix, glmm.genesis$model.matrix, tolerance=1e-6)
    expect_equivalent(nullmod$varComp, glmm.genesis$varComp, tolerance=1e-6)
    expect_equivalent(nullmod$varCompCov, glmm.genesis$varCompCov, tolerance=1e-4)  ## not smaller then 1e-6 because values were not updated after convergence. However, this is not important, because the variance component values are converged at the desired level!
    expect_equivalent(nullmod$family$family, glmm.genesis$family$family)
    expect_equivalent(nullmod$zeroFLAG, glmm.genesis$zeroFLAG)
    expect_true(all(abs(nullmod$cholSigmaInv - glmm.genesis$cholSigmaInv) < 1e-6))
})


test_that("fitNullMod matches fitNullMM - no variance components", {
    n <- 100
    dat <- .testNullInputs(n)
    df <- .asDataFrame(dat)
    
    ## I need to create a junk correlation matrix that will zero out as having variance component zero...  
    varCompJunk <- TRUE
    while(varCompJunk){
	cor.mat <- diag(rnorm(n, sd = 0.01))
	dimnames(cor.mat) <- list(df$scanID, df$scanID)
	lmm.genesis <- GENESIS::fitNullMM(df, "y", covars=c("X1", "X2", "X3"), covMatList=cor.mat, group.var="group", verbose=FALSE)
	if (lmm.genesis$varComp[1] == 0 ) varCompJunk <- FALSE
    }

    nullmod <- fitNullMod(dat$y, dat$X, group.idx=dat$group.idx, verbose=FALSE)

    expect_equivalent(nullmod$fixef, lmm.genesis$fixef, tolerance=1e-6)
    expect_equivalent(nullmod$betaCov, lmm.genesis$betaCov, tolerance=1e-6)
    expect_equivalent(nullmod$resid.marginal, lmm.genesis$resid.marginal, tolerance=1e-6)
    expect_equivalent(nullmod$logLik, lmm.genesis$logLik, tolerance=1e-6)
    expect_equivalent(nullmod$logLikR, lmm.genesis$logLikR, tolerance=1e-6)

    ## currently GENESIS has a mistake, in AUC calculation it uses the number of 
    ## matrices and groups used, but not the actual number of non-zero variance components. 
    ## so this "2" fixes for it. 
    expect_equivalent(nullmod$AIC, lmm.genesis$AIC - 2, tolerance=1e-6)
    expect_equivalent(nullmod$workingY, lmm.genesis$workingY, tolerance=1e-6)
    expect_equivalent(nullmod$model.matrix, lmm.genesis$model.matrix, tolerance=1e-6)
    expect_equivalent(nullmod$varComp, lmm.genesis$varComp[2:3], tolerance=1e-6)
    expect_equivalent(nullmod$varCompCov, lmm.genesis$varCompCov[2:3, 2:3], tolerance=1e-4)
    expect_equivalent(nullmod$family$family, lmm.genesis$family$family)
    expect_true(all(abs(nullmod$cholSigmaInv - lmm.genesis$cholSigmaInv) < 1e-6))
    expect_equivalent(nullmod$RSS, lmm.genesis$RSS, tolerance=1e-6)
})


test_that("nullModelTestPrep matches calculateProjection", {
    n <- 100
    dat <- .testNullInputs(n)
    geno <- .testGenoMatrix(n)
    df <- .asDataFrame(dat)
    
    # basic
    nullmod <- fitNullMod(dat$y, dat$X, verbose=FALSE)
    Xtilde <- calcXtilde(nullmod, geno)
    
    nullmod.orig <- GENESIS::fitNullReg(df, outcome="y", covars=c("X1", "X2", "X3"), verbose=FALSE)
    proj <- GENESIS:::.calculateProjection(nullmod.orig, test="", burden.test="")

    expect_true(all(abs(nullmod$Xtilde - crossprod(proj$Mt, geno)) < 1e-9))
    expect_true(all(abs(nullmod$resid - proj$resid) < 1e-9))
    
    # with covMatList
    nullmod <- fitNullMod(dat$y, dat$X, dat$cor.mat, verbose=FALSE)
    Xtilde <- calcXtilde(nullmod, geno)

    nullmod.orig <- GENESIS::fitNullMM(df, outcome="y", covars=c("X1", "X2", "X3"), covMatList=dat$cor.mat, verbose=FALSE)
    proj <- GENESIS:::.calculateProjection(nullmod.orig, test="", burden.test="")

    expect_true(all(abs(Xtilde - crossprod(proj$Mt, geno)) < 1e-9))
    expect_true(all(abs(nullmod$resid - proj$resid) < 1e-9))
    
    # with group
    nullmod <- fitNullMod(dat$y, dat$X, dat$cor.mat, group.idx=dat$group.idx, verbose=FALSE)
    Xtilde <- calcXtilde(nullmod, geno)

    nullmod.orig <- GENESIS::fitNullMM(df, outcome="y", covars=c("X1", "X2", "X3"), covMatList=dat$cor.mat, group.var="group", verbose=FALSE)
    proj <- GENESIS:::.calculateProjection(nullmod.orig, test="", burden.test="")

    expect_true(all(abs(Xtilde - crossprod(proj$Mt, geno)) < 1e-9))
    expect_true(all(abs(nullmod$resid - proj$resid) < 1e-9))
})


test_that("nullModelTestPrep vs calculateProjection - binary", {
    n <- 100
    geno <- .testGenoMatrix(n)
    
    ## reps <- 0
    ## varCompZero <- TRUE
    ## while(varCompZero & reps < 10){
    
    dat <- .testNullInputs(n, binary=TRUE)
    df <- .asDataFrame(dat)	
    
    ## 	glmm.genesis <- tryCatch({
    ## 		GENESIS::fitNullMM(df, outcome="y", covars=c("X1", "X2", "X3"), covMatList=cor.mat, family="binomial", verbose=FALSE, maxIter=10)
    ## 				}, 
    ## 			warning = function(w){return(list(message = "warning"))},
    ## 			error = function(e){return(list(message = "error"))}
    ## 			)
    ## 	if (!is.null(glmm.genesis$message)) next
    ## 	if (glmm.genesis$varComp[1] != 0 ) varCompZero <- FALSE
    ##         reps <- reps + 1
    ## }
    ## if (varCompZero) stop("could not generate nonzero varComp")

    # basic
    nullmod <- fitNullMod(dat$y, dat$X, family="binomial", verbose=FALSE)
    Xtilde <- calcXtilde(nullmod, geno)
    
    nullmod.orig <- GENESIS::fitNullReg(df, outcome="y", covars=c("X1", "X2", "X3"), family="binomial", verbose=FALSE)
    proj <- GENESIS:::.calculateProjection(nullmod.orig, test="", burden.test="")

    expect_true(all(abs(Xtilde - crossprod(proj$Mt, geno)) < 1e-9))
    expect_true(all(abs(nullmod$resid - proj$resid) < 1e-9))
    
    # with covMatList
    nullmod <- fitNullMod(dat$y, dat$X, dat$cor.mat, family="binomial", verbose=FALSE)
    Xtilde <- calcXtilde(nullmod, geno)

    nullmod.orig <- GENESIS::fitNullMM(df, outcome="y", covars=c("X1", "X2", "X3"), covMatList=dat$cor.mat, family="binomial", verbose=FALSE)
    proj <- GENESIS:::.calculateProjection(nullmod.orig, test="", burden.test="")

    expect_true(all(abs(Xtilde - crossprod(proj$Mt, geno)) < 1e-9))
    expect_true(all(abs(nullmod$resid - proj$resid) < 1e-7))
})
