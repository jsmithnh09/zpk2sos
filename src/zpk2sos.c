#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "zpk2sos.h"

typedef struct {
    double b[3]; // b0, b1, b2
    double a[3]; // a0, a1, a2 (a0 *ought to be* 1.0)
} biquad64_t;

static double get_mag_sq(cplx64_t c) { 
	return c.re * c.re + c.im * c.im; 
}

static double get_dist_sq(cplx64_t c1, cplx64_t c2) {
    double dr = c1.re - c2.re, di = c1.im - c2.im;
    return dr * dr + di * di;
}

/* Comparator for qsort: Sorts poles by distance to unit circle |1 - |p|| */
static int polecmp(const void *a, const void *b) {
    cplx64_t p1 = *(cplx64_t *)a;
    cplx64_t p2 = *(cplx64_t *)b;
    double d1 = fabs(sqrt(get_mag_sq(p1)) - 1.0);
    double d2 = fabs(sqrt(get_mag_sq(p2)) - 1.0);
    if (d1 < d2) return -1;
    if (d1 > d2) return 1;
    return 0;
}


/**
 * core_pairing: Implements nearest-neighbor pairing for padded arrays.
 */
static size_t core_pairing(cplx64_t *z, cplx64_t *p, int n, double k, double *sos_raw) {
    bool *z_used = (bool *)calloc(n, sizeof(bool));
    bool *p_used = (bool *)calloc(n, sizeof(bool));
    int sec_idx = 0;
    int max_stages = soscount(n, n);

    // Pre-sort poles to handle "importance" up front
    qsort(p, n, sizeof(cplx64_t), polecmp);

    for (int i = 0; i < n; i++) {
        if (p_used[i]) continue;

        p_used[i] = true;
        cplx64_t p1 = p[i], p2 = {0, 0};
        bool p2_exists = false;

        // Find conjugate pole
        if (fabs(p1.im) > 1e-12) {
            for (int j = 0; j < n; j++) {
                if (!p_used[j] && fabs(p[j].re - p1.re) < 1e-9 && fabs(p[j].im + p1.im) < 1e-9) {
                    p2 = p[j]; 
		    p_used[j] = true; 
		    p2_exists = true; 
		    break;
                }
            }
        }

        // Find nearest zero to p1
        int z1_idx = -1;
        double min_dist_z = 1e18;
        for (int j = 0; j < n; j++) {
            if (!z_used[j]) {
                double d = get_dist_sq(p1, z[j]);
                if (d < min_dist_z) { 
			min_dist_z = d; 
			z1_idx = j; 
		}
            }
        }

        cplx64_t z1 = z[z1_idx], z2 = {0, 0};
        z_used[z1_idx] = true;
        bool z2_exists = false;

        // Find conjugate zero or second real zero
        if (fabs(z1.im) > 1e-12) {
            for (int j = 0; j < n; j++) {
                if (!z_used[j] && fabs(z[j].re - z1.re) < 1e-9 && fabs(z[j].im + z1.im) < 1e-9) {
                    z2 = z[j]; 
		    z_used[j] = true; 
		    z2_exists = true; 
		    break;
                }
            }
        } else if (p2_exists) {
            int z2_idx = -1; 
	    double min_d2 = 1e18;
            for (int j = 0; j < n; j++) {
                if (!z_used[j] && fabs(z[j].im) < 1e-12) {
                    double d = get_dist_sq(p1, z[j]);
                    if (d < min_d2) { 
			    min_d2 = d; 
			    z2_idx = j; 
		    }
                }
            }
            if (z2_idx != -1) { 
		    z2 = z[z2_idx]; 
		    z_used[z2_idx] = true; 
		    z2_exists = true;
	    }
        }

        // Write to flat double array
        int target = (max_stages - 1 - sec_idx);
	double *s = &sos_raw[target * 6];
        double sk = (sec_idx == max_stages - 1) ? k : 1.0;

        s[0] = sk; // b0
	/* poly expansion. */
        if (z2_exists) {
            s[1] = sk * -(z1.re + z2.re);
            s[2] = sk * (z1.re * z2.re - z1.im * z2.im);
        } else {
            s[1] = sk * -z1.re;
            s[2] = 0.0;
        }

        s[3] = 1.0; // a0
        if (p2_exists) {
            s[4] = -(p1.re + p2.re);
            s[5] = (p1.re * p2.re - p1.im * p2.im);
        } else {
            s[4] = -p1.re;
            s[5] = 0.0;
        }
        sec_idx++;
    }

    free(z_used); 
    free(p_used);
    return (size_t)sec_idx;
}

/**
 * zpk2sos
 * Accepts interleaved doubles [re, im, re, im...]
 * Pads Z and P with zeros at the origin to match lengths.
 */
size_t zpk2sos(const double *z_raw, size_t n_z, const double *p_raw, size_t n_p, double k, double *sos) {
    int max_n = (n_z > n_p) ? n_z : n_p;
    
    cplx64_t *z_pad = (cplx64_t *)calloc(max_n, sizeof(cplx64_t));
    cplx64_t *p_pad = (cplx64_t *)calloc(max_n, sizeof(cplx64_t));

    // Convert interleaved to struct and pad shorter array with (0,0)
    for (size_t i = 0; i < n_z; i++) { 
	    z_pad[i].re = z_raw[2 * i]; 
	    z_pad[i].im = z_raw[2 * i + 1]; 
    }
    for (size_t i = 0; i < n_p; i++) { 
	    p_pad[i].re = p_raw[2 * i]; 
	    p_pad[i].im = p_raw[2 * i + 1]; 
    }

    size_t num_sections = core_pairing(z_pad, p_pad, max_n, k, sos);
    free(z_pad); 
    free(p_pad);
    return num_sections;
}

/**
 * soscount
 * returns the number of biquad sections/stages required for allocation.
 */
size_t soscount(size_t numz, size_t nump) {
	int max_sing = (numz > nump) ? numz : nump;
	return (size_t)((max_sing + 1) / 2); // ceil
}
