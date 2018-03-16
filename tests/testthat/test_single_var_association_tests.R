context("check single variant association tests")

test_that("singleVarTest", {
	n <- 100
	X <- cbind(1, rnorm(n), rbinom(n, size = 1, prob = 0.5))
	y <- X %*% c(1, 0.5, 1) + rnorm(n, sd = c(rep(4, n/2), rep(2, n/2)))
	
	group.idx <- list(G1 = c(1:(n/2)), G2 = c((n/2 + 1):n))
	
	### create a matrix of genetic variants to test.
	geno <- matrix(rbinom(200*n, size = 2, prob = 0.2), nrow = n, ncol = 200)

	nullmod <- fitNullMod(y, X, group.idx = group.idx, verbose=FALSE)
	test.wald <- testGenoSingleVar(nullmod, G = geno, test = "Wald")
	
	# compare to weighted least squares using weighted lm:	
	res.lm <- test.wald
	for (i in 1:ncol(geno)){
		lm.temp <- lm(y ~ -1 + X + geno[,i], weights = c(rep(1/nullmod$varComp[1], n/2), 1/rep(nullmod$varComp[2], n/2)))
		res.lm[i,] <- summary(lm.temp)$coef[4,]
	}
        
        expect_equal(res.lm$Est, test.wald$Est, tolerance = 1e-8)
	expect_equal(res.lm$Wald.Stat, test.wald$Wald.Stat, tolerance = 1e-8)


        # without group
	nullmod <- fitNullMod(y, X, verbose=FALSE)
        test.wald <- testGenoSingleVar(nullmod, G = geno, test = "Wald")
        
        res.lm <- test.wald
	for (i in 1:ncol(geno)){
		lm.temp <- lm(y ~ -1 + X + geno[,i], )
		res.lm[i,] <- summary(lm.temp)$coef[4,]
	}
	
        expect_equal(res.lm$Est, test.wald$Est, tolerance = 1e-8)
	expect_equal(res.lm$Wald.Stat, test.wald$Wald.Stat, tolerance = 1e-8)


        # logistic
        expit <- function(x){exp(x)/(1+exp(x))}
        p <- expit(X %*% c(-1, 0.5, 1))
        D <- rbinom(n, size = 1, prob = p)
		
        ##comparing the Wald test - in genesis computed using the fixed effects of the null model. 
        test.wald <- data.frame(Est = rep(NA, ncol(geno)), Est.SE = NA, Wald.Stat = NA, Wald.pval = NA)
        res.glm <- test.wald	
		
	
	for (i in 1:ncol(geno)){
		nullmod <- fitNullMod(D, cbind(X, geno[,i]), family = "binomial", verbose=FALSE)
		test.wald[i,] <- nullmod$fixef[4,]
                glm.temp <-  glm(D ~ -1 + X + geno[,i], family = "binomial")
                res.glm[i,] <- summary(glm.temp)$coef[4,]
	}

        expect_equal(res.glm$Est, test.wald$Est, tolerance = 1e-8)
	expect_equal(res.glm$Wald.Stat, test.wald$Wald.Stat, tolerance = 1e-8)
	
	## check that we get appropriate error when using the wald test instead of score with binomial outcomes:
	nullmod <- fitNullMod(D, X, family = "binomial", verbose=FALSE)
	expect_message(test.score <- testGenoSingleVar(nullmod, G = geno, test = "Wald"), "Cannot use Wald test")

	expect_equal(colnames(test.score)[1], "Score")
	expect_equal(test.score$Score.pval, test.wald$Wald.pval, tolerance = 0.01)
})


test_that("GxE", {
	n <- 100
	X <- cbind(a=1, b=rnorm(n), c=rbinom(n, size = 1, prob = 0.5))
	y <- X %*% c(1, 0.5, 1) + rnorm(n, sd = c(rep(4, n/2), rep(2, n/2)))
	
	cor.mat <- crossprod(matrix(rnorm(n*n, sd = 0.05),n,n))
        
	### create a matrix of genetic variants to test.
	geno <- matrix(rbinom(200*n, size = 2, prob = 0.2), nrow = n, ncol = 200)

	nullmod <- fitNullMod(y, X, covMatList=cor.mat, verbose=FALSE)
	test.gxe <- testGenoSingleVar(nullmod, G = geno, E = X[,3,drop=FALSE], test = "Wald", GxE.return.cov = TRUE)
        expect_true(all(c("Est.G:c", "SE.G:c") %in% names(test.gxe$res)))
        expect_equal(length(test.gxe$GxEcovMatList), ncol(geno))

        res.lm <- test.gxe$res[,1:4]
        tmp <- data.frame(y, X)
	for (i in 1:ncol(geno)){
		lm.temp <- lm("y ~ b + c + g + c:g", data=cbind(tmp, g=geno[,i]))
		res.lm[i,"Est.G"] <- summary(lm.temp)$coef["g",1]
		res.lm[i,"SE.G"] <- summary(lm.temp)$coef["g",2]
		res.lm[i,"Est.G:c"] <- summary(lm.temp)$coef["c:g",1]
		res.lm[i,"SE.G:c"] <- summary(lm.temp)$coef["c:g",2]
	}
	
        expect_equal(res.lm$Est.G, test.gxe$res$Est.G, tolerance = 1e-8)
        expect_equal(res.lm$SE.G, test.gxe$res$SE.G, tolerance = 1e-8)
        expect_equal(res.lm$`Est.G:c`, test.gxe$res$`Est.G:c`, tolerance = 1e-8)
        expect_equal(res.lm$`SE.G:c`, test.gxe$res$`SE.G:c`, tolerance = 1e-8)

	expect_message(test.gxe <- testGenoSingleVar(nullmod, G = geno, E = X[,3,drop=FALSE], test = "Score"))
})
