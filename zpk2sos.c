
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

typedef struct {
	double re;
	double im;
} cplx64_t;

typedef struct {
	cplx64_t r1;
	cplx64_t r2;
} rootpair;

/**
 * helper functions for complex operations.
 */
static cplx64_t add(cplx64_t a, cplx64_t b)
{
	cplx64_t y = {a.re + b.re, a.im + b.im};
	return y;
}
static cplx64_t sub(cplx64_t a, cplx64_t b)
{
	cplx64_t y = {a.re - b.re, a.im - b.im};
	return y;
}
static cplx64_t mult(cplx64_t a, cplx64_t b)
{
	cplx64_t y = {a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re};
	return y;
}
static double abs64(cplx64_t a)
{
	return sqrt(a.re*a.re + a.im*a.im);
}

static int isreal(cplx64_t r, double tol)
{
	return fabs(r.im) < tol;
}

//validate that a1+b1j and a2-b2j are conjugate; i.e. a1-a2=0 and b1+b2=0.
static int isconj(cplx64_t a, cplx64_t b, double tol)
{
	return fabs(a.re - b.re) < tol && fabs(a.im + b.im) < tol;
}

/**
 * (x - s1)(x - s2) = x^2 - (r1 + r2)x + r1*r2
 * s1 and s2 in this case ought to be conjugates or real. That way the imaginary components
 * get squeezed out.
 */
static void root2quad(cplx64_t s1, cplx64_t s2, double *c0, double *c1, double *c2)
{
	cplx64_t sum = add(s1, s2);
	cplx64_t prod= mult(s1, s2);
	*c0 = 1.0;
	*c1 = -sum.re;
	*c2 = prod.re;
}

/**
 * Imaginary component above is ASSUMED to be non-significant. Otherwise, we have
 * complex coefficients. Alert that we have invalid inputs if this doesn't form a
 * real (non-complex) quadratic.
 */
static int isrealroots(cplx64_t s1, cplx64_t s2)
{
	const double tol = 1e-10;
	cplx64_t sum = add(s1, s2);
	cplx64_t prod = mult(s1, s2);
	if (fabs(sum.im) > tol) return 0;
	if (fabs(prod.im) > tol) return 0;
	return 1;
}

// energy calculation -> \Sigma |x[n]|^2.
double energy(cplx64_t *x, const int N)
{
	double E = 0;
	for (int idx = 0; idx < N; idx++)
	{
		E += abs64(x[idx]) * abs64(x[idx]);
	}
	return E;
}

// distance to unit circle.
static double poledist(cplx64_t p)
{
	return fabs(abs64(p) - 1.0);
}

// distance between two poles.
static double dist(cplx64_t a, cplx64_t b)
{
	double dr = a.re - b.re;
	double di = a.im - b.im;
	return dr*dr + di*di;
}

// follows the same function signature qsort utilizes.
static int polecmp(const void *u, const void *v)
{
	const cplx64_t *a = (const cplx64_t*)u;
	const cplx64_t *b = (const cplx64_t*)v;
	double da = poledist(*a);
	double db = poledist(*b);
	return (da < db) ? -1 : (da > db);
}



// complex conjugate sorting. Doesn't absolute domain cases.
int conjcmp(const void *a, const void *b)
{
    cplx64_t z1 = *(cplx64_t *)a;
    cplx64_t z2 = *(cplx64_t *)b;

    // -1 if ascending, 1 for descending. swap for reversing the sort order.
    if (z1.re != z2.re)
    {
        return (z1.re > z2.re) - (z1.re < z2.re);
    }
    else
    {
        return (z1.im > z2.im) - (z1.im < z2.im);
    }
}

// pre-processing; this will pair roots into conjugate or real pairs first.
static rootpair *pair_roots(const cplx64_t *r, size_t n, size_t *npairs_found)
{
	const double tol = 1e-10;
	int *used = calloc(n, sizeof(int));
	rootpair *pairs = malloc(sizeof(rootpair) * ((n+1)/2));
	size_t npairs = 0;

	for (size_t idx = 0; idx < n; idx++)
	{
		if (used[idx]) continue;

		cplx64_t r1 = r[idx];
		used[idx] = 1;

		/* real root -> pair with itself */
		if (isreal(r1, tol))
		{
			pairs[npairs++] = (rootpair){r1, (cplx64_t){0,0}};
			continue;
		}

		/* find conjugate ahead. */
		int found = 0;
		for (size_t c = idx+1; c < n; c++)
		{
			if (!used[c] && isconj(r1, r[c], tol))
			{
				used[c] = 1;
				pairs[npairs++] = (rootpair){r1, r[c]};
				found = 1;
				break;
			}
		}
		if (!found)
		{
			fprintf(stderr, "ERROR: root (%.6f + %.6fj) has no conjugate.\n",
					r1.re, r1.im);
			free(used);
			free(pairs);
			npairs_found = 0;
			return NULL;
		}
	}
	free(used);
	*npairs_found = npairs;
	return pairs;
}

static int polerootcmp(const void *A, const void *B)
{
	const rootpair *a = (const rootpair*)A;
	const rootpair *b = (const rootpair*)B;

	double da = poledist(a->r1);
	double db = poledist(b->r1);
	return (da < db) ? -1 : (da > db);
}

// function of interest!
size_t zpk2sos(const cplx64_t *z, size_t nz,
		const cplx64_t *p, size_t np,
		double k,
		double *sos)
{
	/* pair the zeros */
	size_t nzpairs = 0;
	rootpair *zpairs = pair_roots(z, nz, &nzpairs);
	if (!nzpairs) return 0;

	/* pair the poles */
	size_t nppairs = 0;
	rootpair *ppairs = pair_roots(p, np, &nppairs);
	if (!nppairs) return 0;

	/* sort pole pairs based on nearest distance. */
	qsort(ppairs, nppairs, sizeof(rootpair), polerootcmp);

	/* number of sections */
	size_t nsec = (nzpairs > nppairs) ? nzpairs : nppairs;
	int *z_used = calloc(nzpairs, sizeof(int));
	for (size_t sing = 0; sing < nsec; sing++)
	{
		/* pick pole pair or pad with real zero (even vs. odd) */
		rootpair pp = (sing < nppairs) ? ppairs[sing] : (rootpair){ {0, 0}, {0,0} };

		/* find nearest zero pair */
		double best = 1e300; // absurd cost.
		int best_idx = -1;

		for (size_t idx = 0; idx < nzpairs; idx++)
		{
			if (!z_used[idx])
			{
				double d = dist(pp.r1, zpairs[idx].r1);
				if (d < best)
				{
					best = d;
					best_idx = (int)idx;
				}
			}
		}
		rootpair zp = (best_idx >= 0) ? zpairs[best_idx] : (rootpair) { {0, 0}, {0, 0} };
		if (best_idx >= 0)
		{
			z_used[best_idx] = 1;
		}

		/* validate we don't have complex coefficients. */
		if (!isrealroots(zp.r1, zp.r2))
		{
			fprintf(stderr, "ERROR: zero pair (%.6f + %.6fj, %.6f + %.6fj) invalid.\n",
					zp.r1.re, zp.r1.im, zp.r2.re, zp.r2.im);
			free(zpairs);
			free(ppairs);
			free(z_used);
			return 0;
		}
		if (!isrealroots(pp.r1, pp.r2))
		{
			fprintf(stderr, "ERROR: pole pair (%.6f + %.6fj, %.6f + %.6fj) invalid.\n",
					pp.r1.re, pp.r1.im, pp.r2.re, pp.r2.im);
			free(zpairs);
			free(ppairs);
			free(z_used);
			return 0;
		}

		/* build the SOS now that we've passed pairing! */
		double b0, b1, b2, a0, a1, a2;
		root2quad(zp.r1, zp.r2, &b0, &b1, &b2);
		root2quad(pp.r1, pp.r2, &a0, &a1, &a2);

		if (sing == 0)
		{
			b0 *= k;
			b1 *= k;
			b2 *= k;
		}

		/* gross pointer math, but it works. */
		double *row = sos + 6*sing;
		row[0]=b0; row[1]=b1; row[2]=b2;
		row[3]=a0; row[4]=a1; row[5]=a2;
	}
	free(zpairs);
	free(ppairs);
	free(z_used);
	return nsec;
}

// for pairing, we check that Nsections = max(ceil(nz/2), ceil(np/2)).
size_t soscount(size_t numzeros, size_t numpoles)
{
	size_t nz = (numzeros + 1) / 2;
	size_t np = (numpoles + 1) / 2;
	return (nz > np) ? nz : np;
}

/*------------------------------- TEST ----------------------------------------------*/
/**
 * poles and zeroes generated by `scipy.signal.ellip(8, 1, 40, 0.2, output='zpk')`.
 * NOTES: SciPy generates the matrix where least sensitive stages come first, most sensitive
 * at the end to prevent numerical precision issues.
 *
 * This is more robust; TODO is to take the result of the zpk2sos and flip the stages. THEN apply the gain.
 */
int main(void)
{
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

	    size_t nz = sizeof(zeros)/sizeof(zeros[0]);
	    size_t np = sizeof(poles)/sizeof(poles[0]);
	    double k = 0.015421175153971328;

	    size_t nsec = soscount(nz, np);
	    printf("Elliptic LPF: %zu zeros, %zu poles → %zu SOS sections\n\n",
        	   nz, np, nsec);

	    double *sos = malloc(nsec * 6 * sizeof(double));
	    zpk2sos(zeros, nz, poles, np, k, sos);

	    for (size_t s = 0; s < nsec; s++) {
	        double *r = sos + 6*s;
        	printf("Section %zu:\n", s+1);
	        printf("  b = [% .6f  % .6f  % .6f]\n", r[0], r[1], r[2]);
        	printf("  a = [% .6f  % .6f  % .6f]\n\n", r[3], r[4], r[5]);
	    }

	    free(sos);
	    return 0;
}
