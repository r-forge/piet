###==========================================================================================================
### Implementation of multiple sample CNV methonds with generalized fused lasso (GFL)
GFL <-
function(Y,Delta,sigma,rho1=1,rho2=2,rho3=0,obj_c=NULL,max_iter=1000,verbose=FALSE) {
    ## Y: #SNP X #SAMPLE
    N <- dim(Y)[1]   ## number of SNPs
    M <- dim(Y)[2]   ## number of sequences
    lambda1 <- rho1 * sigma
    lambda2 <- rho2 * sigma * sqrt(log(N))
    lambda3 <- rho3 * sigma * sqrt(log(N))
        
    obj_c <- ifelse(is.null(obj_c),max(sigma)*0.01,obj_c)
    
    result <- .C("group_fused_lasso",as.double(Y),beta=as.double(Y),as.double(Delta),as.integer(N),as.integer(M),
                                     as.double(lambda1),as.double(lambda2),as.double(lambda3),
                                     obj=as.double(rep(0,max_iter+1)),n_obj=as.integer(0),
                                     as.double(obj_c),as.integer(max_iter),as.integer(verbose))
    
    obj <- result$obj[1:result$n_obj]
    Beta <- matrix(result$beta,N,M)
    return(list(obj=obj,Beta=Beta))
    }
