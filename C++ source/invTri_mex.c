/*
* Description:
* Inverse of triagonal matrix.
*
* based on invChol_mex by Eric Blake
*
* Author(s):
* Gleb Tikhonov
*
* Usage:
* B = invTri_mex(A)
*/

/* If we are on Linux. */
#if !defined(_WIN32) && !defined(_WIN64)
#define dtrtri dtrtri_
#define strtri strtri_
#endif

#include "mex.h"
#include "blas.h"

/* Function declarations */
void dtrtri( char*, char*, mwSize*, double*, mwSize*, mwSignedIndex* );
void strtri( char*, char*, mwSize*, float*, mwSize*, mwSignedIndex* );

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	mwIndex i, j; /* indices for populating lower triangle when finished. */
	mwSize n;   /* Size of the matrix. */
	mwSignedIndex info;     /* flag for success in dpotrf, spotrf, dpotri, spotri */
	char *uplo;     /* upper or lower triangle */
	char *diag = "N";     /* upper or lower triangle */
	mxClassID type;       /* array type */

	/* If we don't pass in the correct number of arguments, throw error. */
	if (nrhs>2 || nrhs==0) {
		mexErrMsgIdAndTxt("mexFunction:invChol_mex:numInputs",
		 "not more than 2 input required: A (square matrix), [UpLo='U' ({'U','L'})]");
	} else if (nrhs==2){
		if(mxIsChar(prhs[1]) != 1){
			mexErrMsgIdAndTxt("MATLAB:invChol_mex:illegaltype", "Second arguments must be either 'U' or 'L' ");
		}
		//uplo = mxCalloc(1,sizeof(char));
		uplo = mxArrayToString(prhs[1]);
		// if(status != 0)
			// mexErrMsgTxt("Cannot allocate memory");
	} else
		uplo = "U";

	type = mxGetClassID(prhs[0]); /* get the array type */

	n = mxGetM(prhs[0]); /* input matrix dimension */
	/* check for symmetric matrix*/
	if (n!=mxGetN(prhs[0])) {
	//	mexErrMsgIdAndTxt("MATLAB:invChol_mex:matchdims",
	//	 "matrix is not symmetric");
	}

	/* check for real matrix */
	if (mxIsComplex(prhs[0])) {
		mexErrMsgIdAndTxt("MATLAB:invChol_mex:iscomplex",
		 "matrix must be real");
	}

	/* create output matrix (fortran modifies the input so we need to copy the input to output) */
	plhs[0]=mxDuplicateArray(prhs[0]);

	/* If we passed in an empty return an empty. */
	if (n==0) {
		return;
	}

	/* double precision */
	if(type==mxDOUBLE_CLASS){
		double *B; /* double pointer to input & output matrices*/
		B = mxGetPr(plhs[0]); /* output matrix pointer */

		dtrtri( uplo, diag, &n, B, &n, &info ); /* Double Cholesky decomposition */
		
		// mexErrMsgIdAndTxt("MATLAB:invChol_mex:dpotri:singularmatrix",
				  // "B %lf", B[3]);
		
		/* check for success */
		if (info<0) {
			mexErrMsgIdAndTxt("MATLAB:invChol_mex:dpotri:illegalvalue",
				  "failed to invert: illegal value");
		}

		if (info>0) {
			mexErrMsgIdAndTxt("MATLAB:invChol_mex:dpotri:singularmatrix",
				  "failed to invert: error in %d", info);
		}

	} else if (type==mxSINGLE_CLASS) {
		float *B; /* float pointer to input and output matrices */
		B = (float*) mxGetData(plhs[0]); /* output matrix pointer */

		strtri( uplo, diag, &n, B, &n, &info ); /* Double Cholesky decomposition */

		/* check for success */
		if (info<0) {
			mexErrMsgIdAndTxt("MATLAB:invChol_mex:spotri:illegalvalue",
				"failed to invert: illegal value");
		}

		if (info>0) {
			mexErrMsgIdAndTxt("MATLAB:invChol_mex:spotri:singularmatrix",
				  "failed to invert: a diagonal element was 0");
		}
		 
	} else {
		mexErrMsgIdAndTxt("MATLAB:invChol_mex:illegaltype",
		 "only single or double matrix inputs are allowed");
	}


	return;
}