###==========================================================================================================
### pre-processing of BAF with or without paired genotype
BAF.transform <- 
function(x,gt=NULL,mBAF.thd=0.97,win.thd=0.8,w=1,k=2,median.adjust=FALSE) {
    require("RANN")
	n <- length(x)   # number of SNPs
	if (median.adjust) {
	    me <- median(x[x>0.25 & x<0.75]) - 0.5
	    x[x>0.25 & x<0.75] <- x[x>0.25 & x<0.75] - me
        }	
	mBAF <- abs(x-0.5) + 0.5
	
	## without paired information
	if (is.null(gt)) {
	    ## mBAF.thd
	    idx <- which(mBAF < mBAF.thd)
	    idx.na <- which(mBAF >= mBAF.thd)
	
	    ## win.thd
	    mBAF1 <- mBAF[idx]
	    n1 <- length(mBAF1)
	    win.sum <- abs(mBAF1-0.5)
	    for (i in 1:w) {
		    win.sum[1:(n1-i)] <- win.sum[1:(n1-i)] + abs(mBAF1[1:(n1-i)]-mBAF1[(i+1):n1])
		    win.sum[(i+1):n1] <- win.sum[(i+1):n1] + abs(mBAF1[(i+1):n1]-mBAF1[1:(n1-i)])
		    }
	    idx1.na <- idx[win.sum >= win.thd]
	    idx.na <- sort(union(idx.na,idx1.na))
	    idx <- sort(setdiff(idx,idx1.na))
		}
	
	## with paired information
	if (!is.null(gt)) {
		idx <- which(gt=="AB")
	    idx.na <- which(gt!="AB")
		}
	
	## impute initial value for missing data
	mBAF1 <- mBAF[idx]
	d.mBAF1 <- mBAF1[-1] - mBAF1[-length(mBAF1)]
	sigma <- sd(d.mBAF1)/sqrt(2)
	
    tmp <- nn2(data=cbind(idx),query=cbind(idx.na),k=k)   ## k-nearest neighbor
    idx.nn <- idx[tmp$nn.idx]
    mBAF[idx.na] <- rowMeans(matrix(mBAF[idx.nn],nrow=length(idx.na))) + rnorm(length(idx.na),0,sigma)
    
	return(list(mBAF=mBAF,idx=idx,idx.na=idx.na))
	}
	