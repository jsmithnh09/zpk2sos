
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdio.h>

// complex conjugate sorting. Doesn't absolute domain cases.
int conjcmp(const void *a, const void *b)
{
    double complex z1 = *(double complex *)a;
    double complex z2 = *(double complex *)b;

    // -1 if ascending, 1 for descending. swap for reversing the sort order.
    if (creal(z1) != creal(z2))
    {
        return (creal(z1) > creal(z2)) - (creal(z1) < creal(z2));
    }
    else
    {
        return (cimag(z1) > cimag(z2)) - (cimag(z1) < cimag(z2));
    }
}

int main(void)
{
    double complex z[4];
    z[0] = 3 + 4*I;
    z[1] = 1 - 1*I;
    z[2] = 2 + 2*I;
    z[3] = 3 - 4*I;
    int n = 4;
    qsort(z, n, sizeof(double complex), conjcmp);

    int idx;
    printf("Testing\n");
    for (idx = 0; idx < n; idx++)
    {
        printf("%lf, %lfi\n", creal(z[idx]), cimag(z[idx]));
    }
    return 0;
}