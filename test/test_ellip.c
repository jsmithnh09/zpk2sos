#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "zpk2sos.h"

/*------------------------------- TEST ---------------------------*/
int main(void)
{
	// this is from the elliptic design; check "gencoeffs.py".
	cplx64_t zeros[] = {
		{0.07501623, 0.99718231}, {0.71403065, 0.70011444},
		{0.78698345, 0.6169741}, {0.80037393, 0.5995011},
        	{0.07501623, -0.99718231}, {0.71403065, -0.70011444},
       		{0.78698345, -0.6169741}, {0.80037393, -0.5995011}
	};
	cplx64_t poles[] = {
		{0.79537728, -0.20999249}, {0.80133831, -0.47682762},
		{0.80453078, -0.56424585}, {0.80647469, -0.58588129},
        	{0.79537728, 0.20999249}, {0.80133831, 0.47682762},
       		{0.80453078, 0.56424585}, {0.80647469, 0.58588129}
	};

	// known SOS.
	double refsos[] = {
		0.01542118, -0.00231368, 0.01542118, 1.000000, -1.59075456, 0.67672186,
		1.00000000, -1.4280613,  1.00000000, 1.000000, -1.60267662, 0.86950766,
		1.00000000, -1.5739669,  1.00000000, 1.000000, -1.60906157, 0.96564316,
		1.00000000, -1.60074787, 1.00000000, 1.000000, -1.61294938, 0.9936583
	};

	    size_t nz = sizeof(zeros)/sizeof(zeros[0]);
	    size_t np = sizeof(poles)/sizeof(poles[0]);
	    double k = 0.015421175153971328;

	    size_t nsec = soscount(nz, np);
	    printf("Elliptic LPF: %zu zeros, %zu poles → %zu SOS sections\n\n",
        	   nz, np, nsec);

	    double *sos = malloc(nsec * 6 * sizeof(double));
	    if (!zpk2sos(zeros, nz, poles, np, k, sos))
	    {
		    printf("Internal error. Exiting.\n");
		    free(sos);
		    return 1;
	    }


	    for (size_t s = 0; s < nsec; s++) 
	    {
	        double *r = sos + 6*s;
        	printf("Section %zu:\n", s+1);
	        printf("  b = [% .6f  % .6f  % .6f]\n", r[0], r[1], r[2]);
        	printf("  a = [% .6f  % .6f  % .6f]\n\n", r[3], r[4], r[5]);
	    }
	    for (size_t s = 0; s < nsec; s++)
	    {
		double *r = sos + 6*s;
		double *v = refsos + 6*s;
		for (size_t c = 0; c < 6; c++)
		{
			if (fabs(r[c] - v[c]) > 1e-6)
			{
				printf("FAIL precision on section %zu.\n", s);
				free(sos);
				return 1;
			}
		}
	    }
	    printf("PASS\n");
	    free(sos);
	    return 0;
}

