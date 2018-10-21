

/* If we are on Linux. */
#if !defined(_WIN32) && !defined(_WIN64)
#define dpotrf dpotrf_
#define spotrf spotrf_
#define dpotri dpotri_
#define spotri spotri_
#endif

#include "mex.h"
#include "blas.h"

/* Function declarations */
void dpotrf( char*, mwSize*, double*, mwSize*, mwSignedIndex* );
void dpotri( char*, mwSize*, double*, mwSize*, mwSignedIndex* );
void spotrf( char*, mwSize*, float*, mwSize*, mwSignedIndex* );
void spotri( char*, mwSize*, float*, mwSize*, mwSignedIndex* );

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    mwIndex i, j, k; /* indices for populating lower triangle when finished. */
	 mwSize *dimVec; // dimensions of the array
    mwSize n, K;   /* Size of the matrix. */
    mwSignedIndex info;     /* flag for success in dpotrf, spotrf, dpotri, spotri */
    char *uplo = "L";     /* upper or lower triangle */
    mxClassID type;       /* array type */
    
    /* If we don't pass in the correct number of arguments, throw error. */
    if (nrhs!=1) {
        mexErrMsgIdAndTxt("mexFunction:invCholx_mex:numInputs",
                "1 input required: A (square matrix)");
    }
    
    type = mxGetClassID(prhs[0]); /* get the array type */
    
	 dimVec = mxGetDimensions(prhs[0]); /* input array dimension */
	 n = dimVec[0];
	 K = dimVec[2];
	 char buffer [50];
	 // sprintf(buffer, "%d %d %d", dimVec[0], dimVec[1], dimVec[2]);
	 // mexErrMsgIdAndTxt("MATLAB:invCholx_mex:matchdims",
                // buffer);
    // n = mxGetM(prhs[0]); /* input matrix dimension */
    // /* check for symmetric matrix*/
    if (n!=dimVec[1]) {
        mexErrMsgIdAndTxt("MATLAB:invCholx_mex:matchdims",
                "matrix is not symmetric");
    }
	
	/* check for real matrix */
	if (mxIsComplex(prhs[0])) {
		mexErrMsgIdAndTxt("MATLAB:invCholx_mex:iscomplex",
                "matrix must be real");
	}
	
    /* create output matrix (fortran modifies the input so we need to copy the input to output) */
    plhs[0]=mxDuplicateArray(prhs[0]);
    
    /* If we passed in an empty return an empty. */
    if (n==0) {
        return;
    }
    
    /* double precision */
    if(type==mxDOUBLE_CLASS)
    {
		 for(k=0; k<K; ++k){
        double *B; /* double pointer to input & output matrices*/
        B = mxGetPr(plhs[0]) + k*n*n; /* output matrix pointer */
              
        dpotrf( uplo, &n, B, &n, &info ); /* Double Cholesky decomposition */
        
        /* check for success */
        if (info<0) {
            mexErrMsgIdAndTxt("MATLAB:invChol_mex:dpotrf:illegalvalue",
                    "cholesky decomposition failed: illegal value ");
        }
        if (info>0) {
            mexErrMsgIdAndTxt("MATLAB:invChol_mex:dpotrf:notposdef",
                    "cholesky decomposition failed: matrix is not positive definite");
        }
		  for (i=0; i<n; i++) {
            for (j=i+1; j<n; j++) {
                B[j*n+i]=0;
            }
        }
		 }
        
    } else if (type==mxSINGLE_CLASS) {
		 for(k=0; k<K; ++k){
        float *B; /* float pointer to input and output matrices */
        B = (float*) mxGetData(plhs[0]) + sizeof(float)*k*n*n; /* output matrix pointer */
        
        spotrf( uplo, &n, B, &n, &info ); /* Double Cholesky decomposition */
        
        /* check for success */
        if (info<0) {
            mexErrMsgIdAndTxt("MATLAB:invChol_mex:spotrf:illegalvalue",
                    "cholesky decomposition failed: illegal value ");
        }
        if (info>0) {
            mexErrMsgIdAndTxt("MATLAB:invChol_mex:spotrf:notposdef",
                    "cholesky decomposition failed: matrix is not positive definite");
        }
        for (i=0; i<n; i++) {
            for (j=i+1; j<n; j++) {
                B[j*n+i]=0;
            }
        }
      
		 }
        
    } else {
        mexErrMsgIdAndTxt("MATLAB:invChol_mex:illegaltype",
                "only single or double matrix inputs are allowed");
    }
    
    
    return;
}