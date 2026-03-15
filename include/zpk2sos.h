/**
 * Zero Pole Gain to Second Order Section Conversion.
 * Author: Jordan R. Smith, 2026.
 */

#ifndef _ZPK2SOS_H_
#define _ZPK2SOS_H_

#ifdef __cplusplus
extern "C" {
#endif


// expose the functions for use.
#ifdef _WIN32
#define EXPORT __declspec(dllexport)
#else
#define EXPORT __attribute__((visibility("default")))
#endif

#include <stdlib.h>

/**
 * struct cplx64_t - simplified complex number interface.
 * The windows API is needlessly confusing, and although Linux has
 * it correct, I'd rather just define the basic operations upfront.
 *
 * @re: real component of a complex number.
 * @im: imaginary component of a complex number.
 */
typedef struct { double re; double im; } cplx64_t;

/**
 * zpk2sos() -- group nearest poles and zeros together, matching the "nearest"
 * method approach used in SciPy.
 * @z:   an array of complex zeroes from filter design.
 * @numz:   the number of zeroes in the z-array.
 * @p:   an array of complex poles for pairing.
 * @nump:   the number of poles in the p-array.
 * @k:   the real gain associated with the system.
 * @sos:   the pre-allocated SOS matrix of size (numsections, 6).
 *
 * Return: non-zero integer indicating success; zero indicating an error.
 */
size_t _zpk2sos(const cplx64_t *z, size_t numz, const cplx64_t *p, size_t nump, double k, double *sos);

/**
 * zpk2sos() -- same function as above, but instead expects the zeros to be
 * real and imaginary interlaced doubles.
 * @z:	an array of real/imaginary double zeroes from filter design, (2x numz)
 * @p:  an array of real/imaginary double poles from filter design, (2x nump)
 * @k:	the real gain associated with the system.
 * @sos:  the pre-allocated SOS matrix of size (numsections, 6).
 *
 * Return: non-zero integer indicating success; zero indicating error.
 */
EXPORT size_t zpk2sos(const double *z, size_t numz, const double *p, size_t nump, double k, double *sos);

/**
 * soscount() -- returns the total number of required biquads.
 * @numzeros:   indicates the number of zero singularities.
 * @numpoles:   indicates the number of pole singularities.
 *
 * Return: the number of biquads required in the SOS matrix.
 */
EXPORT size_t soscount(size_t numzeros, size_t numpoles);


#ifdef __cplusplus
}
#endif

#endif /* _ZPK2SOS_H_ */
