#include <math.h>
#include <string.h>
#include "mex.h"

/***************************************************************************
 * WARNING: NO ERROR CHECKING WHATSOEVER - CALL IT RIGHT!!! (see genvar.m) *
 ***************************************************************************
 *
 *     X = genvar_mex(A,Z);
 *
 */

void mexFunction(const int nlhs, mxArray *plhs[], const int nrhs, const mxArray const *prhs[])
{
	const double* const a  = mxGetPr(prhs[0]);
	const double* const z  = mxGetPr(prhs[1]);
	
	const size_t n = mxGetM(prhs[1]); 
	const size_t m = mxGetN(prhs[1]);
	const mwSize d = mxGetNumberOfDimensions(prhs[0]);
	const mwSize p = (d < 3 ? 1 : mxGetDimensions(prhs[0])[2]);
	
	double* const x = mxGetPr(plhs[0] = mxCreateDoubleMatrix(n,m,mxREAL));/* allocate output variable x */
	memcpy(x,z,n*m*sizeof(double));                                       /* and copy it from z         */

	const size_t nn = n*n; 
	size_t t;
	double* xt = x+n*p;
	for (t=p; t<m; ++t) {
		size_t k;
		const double* xtk = xt-n;
		const double* ak = a;
		for (k=0; k<p; ++k) {
			size_t i;
			double* xti = xt;
			const double* aki = ak;
			for (i=0; i<n; ++i) {
				size_t j;
				double w = 0.0;
				for (j=0; j<n; ++j) w += *(aki+n*j) * *(xtk+j);
				*xti++ += w;
				++aki;
			}
			xtk -= n;
			ak += nn;
		}
		xt += n;
	}
}
