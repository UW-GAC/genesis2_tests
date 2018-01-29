
## takes a null model and prepare specific arguments to streamline the testing
nullModelTestPrep <- function(nullmod, G){
    
    Y <- nullmod$workingY
    W <- nullmod$model.matrix
   
    if (nullmod$family$mixedmodel){  ## n by n cholSigmaInv
        C <- nullmod$cholSigmaInv
        CW <- crossprod(C, W)
        MtG <- crossprod(C,G) - tcrossprod(tcrossprod(CW, chol2inv(chol(crossprod(CW)))), crossprod(G,tcrossprod(C,t(CW))))
        rm(G)
        Ytilde <- crossprod(C,Y) - crossprod(tcrossprod(tcrossprod(C, tcrossprod(chol2inv(chol(crossprod(CW))), CW)), CW),Y)
        resid <- crossprod(t(C),Ytilde) - crossprod(t(tcrossprod(tcrossprod(C, tcrossprod(chol2inv(chol(crossprod(CW))), CW)), CW)),Ytilde)
        
        # Projection matrix P = Mt %*% M
        # Ytilde = phenotype adjusted for the covariates/correlation structure
        #Mt <- C - tcrossprod(tcrossprod(C, tcrossprod(chol2inv(chol(crossprod(CW))), CW)), CW)
        #rm(C)
        #Ytilde <- crossprod(Mt, Y)
        #resid <- as.vector(Mt %*% crossprod(Mt, Y))
    }
    
    if (!nullmod$family$mixedmodel & (nullmod$family$family != "gaussian")){
        ## above we use cholSigmaInv, here we use Sigma, but the math is the same.
        ## why is this case different from the gaussian case below?
        sigma <- sqrt(nullmod$varComp)
        C <- Diagonal(x=sigma)
        CW <- W * sigma
        MtG <- crossprod(C,G) - tcrossprod(tcrossprod(CW, chol2inv(chol(crossprod(CW)))), crossprod(G,tcrossprod(C,t(CW))))
        rm(G)
        Ytilde <- crossprod(C,Y) - crossprod(tcrossprod(tcrossprod(C, tcrossprod(chol2inv(chol(crossprod(CW))), CW)), CW),Y)
        
        # Projection matrix P = Mt %*% M
        # Ytilde = phenotype adjusted for the covariates/correlation structure
        #Mt <- C - tcrossprod(tcrossprod(C, tcrossprod(chol2inv(chol(crossprod(CW))), CW)), CW)
        #rm(C)
        #Ytilde <- crossprod(Mt, Y)
        resid <- nullmod$resid.marginal
    }

    
    if (!nullmod$family$mixedmodel & (nullmod$family$family == "gaussian")){  ## a diagonal or scalar cholSigmaInv
        
        if (nullmod$hetResid)	{  ## cholSigmaInv is diagonal
            C <- diag(nullmod$cholSigmaInv)
        } else { ## family is "gaussian", cholSigmaInv is a scalar.
            C <- nullmod$cholSigmaInv        
        }	
        
        CW <- W * C      ## this is equal to crossprod(diag(C), W) when C is a vector  
        MtG <- G*C - tcrossprod(tcrossprod(CW, chol2inv(chol(crossprod(CW)))), crossprod(G,CW*C))
        rm(G)
        Ytilde <- C*Y - crossprod(tcrossprod(t(tcrossprod(chol2inv(chol(crossprod(CW))), CW))*C, CW) ,Y)
        
        #Mt <- -tcrossprod(t(tcrossprod(chol2inv(chol(crossprod(CW))), CW))*C, CW)
        #diag(Mt) <- diag(Mt) + C
        #rm(C)
        #Ytilde <- crossprod(Mt, Y)
        
        ## prepare resids for testing
        if (nullmod$hetResid){
            resid <- as.vector(C*Ytilde - crossprod(t(tcrossprod(t(tcrossprod(chol2inv(chol(crossprod(CW))), CW))*C, CW)),as.vector(Ytilde)))
            #resid <- as.vector(Mt %*% crossprod(Mt, Y))
        } else{
            resid <- nullmod$resid.marginal/nullmod$varComp
        }
        
    }	
    
    
    
    sY2 <- sum(Ytilde^2)

    out <- list(Xtilde = MtG, Ytilde = Ytilde, sY2 = sY2, k = ncol(W), resid = resid, family = nullmod$family$family)
    return(out)
}


## idx.exclude are indices of individuals that should be excluded (e.g. because of missing genotypes)
nullModelSubset <- function(nullmod, idx.exclude){
    for (v in c("sample.id", "fitted.values", "resid.marginal", "resid.condition", "workingY", "outcome")) {
        nullmod[[v]] <- nullmod[[v]][-idx.exclude]
    }
    nullmod$model.matrix <- nullmod$model.matrix[-idx.exclude,]
    
    if (nullmod$family$mixedmodel){  ## n by n cholSigmaInv {
        nullmod$cholSigmaInv <- subsetCholSigmaInv(nullmod$cholSigmaInv, idx.exclude)
    }
    
    if (!nullmod$family$mixedmodel & (nullmod$family$family == "gaussian")){  ## a diagonal or scalar cholSigmaInv
        
        if (nullmod$hetResid)	{  ## cholSigmaInv is diagonal
            nullmod$cholSigmaInv <- nullmod$cholSigmaInv[-idx.exclude, -idx.exclude]
        }   
    }
}


# this is a fancy way of getting the inverse of the subset without having to get the original matrix
# cholesky decomposition of sigma inverse (inverse phenotype covariance matrix)
subsetCholSigmaInv <- function(cholSigmaInv, chol.idx) {
    if(length(chol.idx) > 0){
        # subset cholSigmaInv
        SigmaInv <- tcrossprod(cholSigmaInv)
        for(i in sort(chol.idx, decreasing=TRUE)){
            SigmaInv <- SigmaInv[-i,-i] - tcrossprod(SigmaInv[-i,i])/SigmaInv[i,i]
        }
        cholSigmaInv <- t(chol(SigmaInv))
    }
    
    cholSigmaInv
}
