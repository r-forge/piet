###==========================================================================================================
### Implementation of fused lasso for single sample
FL <- 
function(y,sigma,rho1=1,rho2=2,obj_c=1e-4,max_iter=1000) {
    N <- length(y)   ## number of SNPs
    lambda1 <- rho1 * sigma
    lambda2 <- rho2 * sigma * sqrt(log(N))

    ##beta <- y   ## latent copy number for each sample each SNP, initialize

    result <- .C("triDiagonal",as.double(y),beta=as.double(y),as.integer(N),
                               as.double(lambda1),as.double(lambda2),
                               obj=as.double(rep(0,max_iter)),n_obj=as.integer(0),
                               as.double(obj_c),as.integer(max_iter) ) 

    return(list(obj=result$obj[1:result$n_obj],beta=result$beta))
    }