#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>

#define eps 1e-10

#define DEBUG 1

void triDiagonal(const double *y, double *x, const int *n,
				 const double *l1, const double *l2, double *obj, int *n_obj,
				 const double *obj_c, const int *max_iter);

void update_individual(const double *y, double *x, int *n,const double *gfl_coef, 
					   const double *l1, const double *l2, const double *l3);

void update_individual1(const double *y, double *x, int *n, const double *delta, const double *gfl_coef, 
					   const double *l1, const double *l2, const double *l3);

void segment(const double *x, const int *n, const double *min_step, const int *min_len, int *segments, int *s);


void group_fused_lasso(const double *signal, double *beta, const double *idxNA,  const int *n, const int *m, 
					   const double *lambda1, const double *lambda2, const double *lambda3,
					   double *obj, int *n_obj,const double *obj_c, const int *max_iter, const int *verbose);



void group_fused_lasso(const double *signal, double *beta, const double *idxNA,  const int *n, const int *m, 
					   const double *lambda1, const double *lambda2, const double *lambda3,
					   double *obj, int *n_obj,const double *obj_c, const int *max_iter, const int *verbose)
{
	// Note: R adopt column major rule
	// signal: original signal for each individual (in vector n*m)
	// beta: the fitted piece-wise constant function for each individual (in vector n*m)
	// idxNA: idxNA=1, observation exists at that locus
	//        idxNA=0. missing data at that locus (in vector n*m)
	// n: the number of SNPs
	// m: the number of samples	
	// coef3 is group-fused-lasso penaly coefficient for each SNP across all individuals
	// lambda1, lambda2 and lambda3 are tuning parameters
	// obj: the vector of the values of objective function at each iteration
	// n_obj: the length of the vector obj
	// obj_c: convergence criterion defined in terms of the improvement of objective function
	// max_iter: maximum allowed iteration
	// verbose: whether display some imtermediate diagnosis information. 1=Yes, 0=No


	double loss,penalty1,penalty2,penalty3;
	double loss1,p1,p2,p3;

	double *dA = calloc(*n,sizeof(double));
	double *B = calloc(*n,sizeof(double));
	double *dA1 = calloc(*n-1,sizeof(double));
	double *coef3 = calloc(*n-1,sizeof(double));
	double *temp_coef = calloc(*n-1,sizeof(double));
	double temp;
	
	// data are transposed: #Sample X #SNP
	double **Y = (double **) calloc(*m,sizeof(double *));
	double **X = (double **) calloc(*m,sizeof(double *));
	double **Delta = (double **) calloc(*m,sizeof(double *));	
	
	for (int i=0; i<*m; i++) {
		Y[i] = (double *) calloc(*n,sizeof(double));
		X[i] = (double *) calloc(*n,sizeof(double));
		Delta[i] = (double *) calloc(*n,sizeof(double));
	}
		
	// transform the data vector input by R into matrix
	for (int i=0; i<*m; i++) {
		for (int j=0; j<*n; j++) {
			Y[i][j] = signal[(*n)*i+j];
			X[i][j] = beta[(*n)*i+j];
		    Delta[i][j] = idxNA[(*n)*i+j];
		}
	}
	
	// initialize object function
	loss = 0;
	penalty1 = 0;
	penalty2 = 0;
	penalty3 = 0;
	
	// compute loss, penalty1
	for (int i=0; i<*m; i++) {
		loss1 = 0;
		p1 = 0;
		for (int j=0; j<*n; j++) {
			temp = Y[i][j]*Delta[i][j] - X[i][j]*Delta[i][j];
			loss1 += temp*temp;
		    p1 += sqrt(X[i][j]*X[i][j] + eps);
		}
		loss += 0.5 * loss1;
		penalty1 += lambda1[i] * p1;
	}
	// compute penalty2
	for (int i=0; i<*m; i++) {
		p2 = 0;
		for (int j=0; j<*n-1; j++) {
			temp = X[i][j+1] - X[i][j];
		    p2 += sqrt(temp*temp + eps);
		}
		penalty2 += lambda2[i] * p2;
	}
	// compute penalty3
	for (int j=0; j<*n-1; j++) {
		p3 = 0;
		for (int i=0; i<*m; i++) {
			temp = lambda3[i] * (X[i][j+1] - X[i][j]);
		    p3 += temp*temp;
		}
		coef3[j] = sqrt(p3 + eps);
		penalty3 += coef3[j];
	}
	obj[0] = loss + penalty1 + penalty2 + penalty3;
	
	if (*verbose==1) printf("%4d, %15.6f\n",1,obj[0]);
	
	
	for (int k=0; k<*max_iter; k++) {   // kth iteration
		
		for (int i=0; i<*m; i++) {   // ith sample
			
			// update each sample
			// Update the linear system A*beta=B to solve the quadratic optimization problem
			// Update matrix A
			// Update diagonal entries	
			
			for (int j=0; j<*n-1; j++) {
				temp = X[i][j+1] - X[i][j];
				temp_coef[j] = lambda2[i]/sqrt(temp*temp+eps) + lambda3[i]*lambda3[i]/coef3[j];
			}
			
			for (int j=0; j<*n; j++)    dA[j] = Delta[i][j] + lambda1[i]/sqrt(X[i][j]*X[i][j]+eps);
			for (int j=0; j<*n-1; j++)  dA[j] += temp_coef[j];
			for (int j=1; j<*n; j++)    dA[j] += temp_coef[j-1];
			
			// Update 1-order off-diagonal entries
			
			for (int j=0; j<*n-1; j++) dA1[j] = -temp_coef[j];
			
			// Updata vector B
			
			for (int j=0; j<*n; j++) B[j] = Y[i][j] * Delta[i][j];
			
			// Apply "chasing" algorithm to solve this linear system
			// Forward sweeping 1-order off-diagonal entries
			
			for (int j=0; j<*n-1; j++) {
				dA[j+1] -= dA1[j]*dA1[j]/dA[j];
				B[j+1] -= B[j]*dA1[j]/dA[j];
			}
			
			// Backward sweeping 1-order off-diagonal entries
			
			for (int j=0; j<*n-1; j++) B[*n-2-j] -= B[*n-1-j]*dA1[*n-2-j]/dA[*n-1-j];
			
			// Update beta
			
			for (int j=0; j<*n; j++) X[i][j] = B[j]/dA[j];
			
		} // for i
		
	
		// compute object function
		loss = 0;
		penalty1 = 0;
		penalty2 = 0;
		penalty3 = 0;
		
		// compute loss, penalty1
		for (int i=0; i<*m; i++) {
			loss1 = 0;
			p1 = 0;
			for (int j=0; j<*n; j++) {
				temp = Y[i][j]*Delta[i][j] - X[i][j]*Delta[i][j];
				loss1 += temp*temp;
				p1 += sqrt(X[i][j]*X[i][j] + eps);
			}
			loss += 0.5 * loss1;
			penalty1 += lambda1[i] * p1;
		}
		// compute penalty2
		for (int i=0; i<*m; i++) {
			p2 = 0;
			for (int j=0; j<*n-1; j++) {
				temp = X[i][j+1] - X[i][j];
				p2 += sqrt(temp*temp + eps);
			}
			penalty2 += lambda2[i] * p2;
		}
		// compute penalty3
		for (int j=0; j<*n-1; j++) {
			p3 = 0;
			for (int i=0; i<*m; i++) {
				temp = lambda3[i] * (X[i][j+1] - X[i][j]);
				p3 += temp*temp;
			}
			coef3[j] = sqrt(p3 + eps);
			penalty3 += coef3[j];
		}
		obj[k+1] = loss + penalty1 + penalty2 + penalty3;
		
		if (*verbose==1) printf("%4d, %15.6f\n",k+2,obj[k+1]);
		
		if ( fabs(obj[k+1]-obj[k]) < *obj_c ) {
			n_obj[0] = k+2;
			break;
		}
		
	} // for k
	
	// transform the data matrix into vecter output to R
	for (int i=0; i<*m; i++) {
		for (int j=0; j<*n; j++) {
		    beta[(*n)*i+j] = X[i][j];
		}
	}
}


void update_individual1(const double *y, double *x, int *n, const double *delta, const double *gfl_coef, 
						const double *l1, const double *l2, const double *l3)
{
	// y is original signal for each individual
	// x is the beta for each individual
	// n the number of SNPs
	// delta: delta=1, observation exists for that locus
	//        delta=0. missing data for that locus
	// gfl_coef is group-fused-lasso penaly coefficient for each SNP across all individuals
	// l1, l2 and l3 are tuning parameters
    
	double dA[*n],B[*n];
	double dA1[*n-1],temp_coef[*n-1];
	double temp;
	
	// Update the linear system A*beta=B to solve the quadratic optimization problem
	// Update matrix A
	// Update diagonal entries	
	
	for (int i=0; i<*n-1; i++) {
		temp = x[i+1] - x[i];
		temp_coef[i] = *l2/sqrt(temp*temp+eps) + (*l3)*(*l3)/gfl_coef[i];
	}
	
	for (int i=0; i<*n; i++)    dA[i] = delta[i] + (*l1)/sqrt(x[i]*x[i]+eps);
	for (int i=0; i<*n-1; i++)  dA[i] = dA[i] + temp_coef[i];
	for (int i=1; i<*n; i++)    dA[i] = dA[i] + temp_coef[i-1];
	
	// Update 1-order off-diagonal entries
	
    for (int i=0; i<*n-1; i++) dA1[i] = -temp_coef[i];
	
	// Updata vector B
	
	for (int i=0; i<*n; i++) B[i] = y[i];
	
	// Apply "chasing" algorithm to solve this linear system
	// Forward sweeping 1-order off-diagonal entries
	
	for (int i=0; i<*n-1; i++) {
		dA[i+1] = dA[i+1] - dA1[i]*dA1[i]/dA[i];
		B[i+1] = B[i+1] - B[i]*dA1[i]/dA[i];
	}
	
	// Backward sweeping 1-order off-diagonal entries
	
	for (int i=0; i<*n-1; i++) B[*n-2-i] = B[*n-2-i] - B[*n-1-i]*dA1[*n-2-i]/dA[*n-1-i];
	
	// Update beta
	
	for (int i=0; i<*n; i++) x[i] = B[i]/dA[i];
}


void triDiagonal(const double *y, double *x, const int *n,
				 const double *l1, const double *l2, double *obj, int *n_obj,
				 const double *obj_c, const int *max_iter)
{
	// y is original signal for each individual
	// x is the beta for each individual
	// n the number of SNPs
	// l1 and l2 are tuning parameters
	
	double *dA = (double *) calloc(*n,sizeof(double));
	double *B = (double *) calloc(*n,sizeof(double));
	double *dA1 = (double *) calloc(*n-1,sizeof(double));
	double *temp_coef = (double *) calloc(*n-1,sizeof(double));
	double temp;
	double loss,penalty1,penalty2;	
	
	loss = 0;
	penalty1 = 0;
	for (int i=0; i<*n; i++) {
		temp = y[i] - x[i];
		loss += temp*temp;
		penalty1 += sqrt(x[i]*x[i] + eps);
	}
	penalty2 = 0;
	for (int i=0; i<*n-1; i++) {
		temp = x[i] - x[i+1];
		penalty2 += sqrt(temp*temp + eps);
	}
	obj[0] = 0.5*loss + (*l1)*penalty1 + (*l2)*penalty2;

		
	for (int k=0; k<*max_iter; k++) {

		//for (int i=0; i<*n; i++) old_x[i] = x[i];
		
		// Update the linear system A*beta=B to solve the quadratic optimization problem
		// Update matrix A
		// Update diagonal entries	
	
		for (int i=0; i<*n-1; i++) {
			temp = x[i+1] - x[i];
			temp_coef[i] = (*l2)/sqrt(temp*temp+eps);
		}
	
		for (int i=0; i<*n; i++)    dA[i] = 1 + *l1/sqrt(x[i]*x[i]+eps);
		for (int i=0; i<*n-1; i++)  dA[i] = dA[i] + temp_coef[i];
		for (int i=1; i<*n; i++)    dA[i] = dA[i] + temp_coef[i-1];
	
		// Update 1-order off-diagonal entries
	
	    for (int i=0; i<*n-1; i++) dA1[i] = -temp_coef[i];
	
		// Updata vector B
	
		for (int i=0; i<*n; i++) B[i] = y[i];
	
		// Apply "chasing" algorithm to solve this linear system
		// Forward sweeping 1-order off-diagonal entries
	
		for (int i=0; i<*n-1; i++) {
			dA[i+1] = dA[i+1] - dA1[i]*dA1[i]/dA[i];
			B[i+1] = B[i+1] - B[i]*dA1[i]/dA[i];
		}
	
		// Backward sweeping 1-order off-diagonal entries
	
		for (int i=0; i<*n-1; i++) B[*n-2-i] = B[*n-2-i] - B[*n-1-i]*dA1[*n-2-i]/dA[*n-1-i];
	
		// Update beta
	
		for (int i=0; i<*n; i++) x[i] = B[i]/dA[i];
		

		loss = 0;
		penalty1 = 0;
		for (int i=0; i<*n; i++) {
			temp = y[i] - x[i];
			loss += temp*temp;
			penalty1 += sqrt(x[i]*x[i] + eps);
		}
		penalty2 = 0;
		for (int i=0; i<*n-1; i++) {
			temp = x[i] - x[i+1];
			penalty2 += sqrt(temp*temp + eps);
		}
		obj[k+1] = 0.5*loss + (*l1)*penalty1 + (*l2)*penalty2;
		
		if ( fabs(obj[k+1]-obj[k]) < *obj_c ) {
			n_obj[0] = k+2;
			break;
		}
    }
}




void update_individual(const double *y, double *x, int *n,const double *gfl_coef, 
					   const double *l1, const double *l2, const double *l3)
{
// y is original signal for each individual
// x is the beta for each individual
// n the number of SNPs
// gfl_coef is group-fused-lasso penaly coefficient for each SNP across all individuals
// l1, l2 and l3 are tuning parameters
    
	double dA[*n],B[*n];
	double dA1[*n-1],temp_coef[*n-1];
	double temp;
	
// Update the linear system A*beta=B to solve the quadratic optimization problem
// Update matrix A
// Update diagonal entries	

	for (int i=0; i<*n-1; i++) {
		temp = x[i+1] - x[i];
		temp_coef[i] = *l2/sqrt(temp*temp+eps) + (*l3)*(*l3)/gfl_coef[i];
	}
	
	for (int i=0; i<*n; i++)    dA[i] = 1 + *l1/sqrt(x[i]*x[i]+eps);
	for (int i=0; i<*n-1; i++)  dA[i] = dA[i] + temp_coef[i];
	for (int i=1; i<*n; i++)    dA[i] = dA[i] + temp_coef[i-1];
	
// Update 1-order off-diagonal entries

    for (int i=0; i<*n-1; i++) dA1[i] = -temp_coef[i];

// Updata vector B

	for (int i=0; i<*n; i++) B[i] = y[i];

// Apply "chasing" algorithm to solve this linear system
// Forward sweeping 1-order off-diagonal entries

	for (int i=0; i<*n-1; i++) {
		dA[i+1] = dA[i+1] - dA1[i]*dA1[i]/dA[i];
		B[i+1] = B[i+1] - B[i]*dA1[i]/dA[i];
	}
	
// Backward sweeping 1-order off-diagonal entries

	for (int i=0; i<*n-1; i++) B[*n-2-i] = B[*n-2-i] - B[*n-1-i]*dA1[*n-2-i]/dA[*n-1-i];
	
// Update beta

	for (int i=0; i<*n; i++) x[i] = B[i]/dA[i];
}




void segment(const double *x, const int *n, const double *min_step, const int *min_len, int *segments, int *s)
{
// x: a squence to be segmented
// n: the length of n
// min_step: minimum difference between pieces
// min_len: minimum length of a piece
// return: 
// segments: a list of segments with start and end idx
// s: the number of segments
	
	int left,right,ns;

	ns = 0;
	left = 0;
	right = left + 1;
	
	while (left <= *n - *min_len) {
		
		if ( fabs(x[right]-x[right-1]) <= *min_step ) {
			if ( right < *n - *min_len ) {
				right++;
			    continue;
			}
			if ( right >= *n - *min_len ) {
				ns++;
				segments[2*(ns-1)] = left+1;
				segments[2*(ns)-1] = *n;
				break;
			}
		}
		
		if ( fabs(x[right]-x[right-1]) > *min_step ) {
			if ( right < *n - *min_len ) {
				if ( right-left >= *min_len) {
					ns++;
					segments[2*(ns-1)] = left + 1;
					segments[2*(ns)-1] = right;
					left = right;
					right = left + 1;
					continue;
				}
				if ( right-left < *min_len ) {
					right++;
					continue;
				}
			}
			
			if (right >= *n - *min_len) {
				ns++;
				segments[2*(ns-1)] = left+1;
				segments[2*(ns)-1] = *n;
				break;
			}
		}
	}
	s[0] = ns;
}
