#include <math.h>
#include "groups.hpp"
#include "nrutil.h"
#define EPS 1.0e-6
#define JMAX 14
#define JMAXP (JMAX+1)
#define K 5

/**
 * @file qromo.c
 * @brief Numerical integration using Romberg's method.
 *
 * This file contains the implementation of the qromo function, which performs
 * numerical integration of a given function over a specified interval using
 * Romberg's method.
 *
 * @param func Pointer to the function to be integrated. The function should take
 *             a single float argument and return a float.
 * @param a    The lower limit of integration.
 * @param b    The upper limit of integration.
 * @param choose Pointer to a function that selects the appropriate integration
 *               method. This function should take a pointer to the function to
 *               be integrated, the lower and upper limits of integration, and
 *               an integer parameter, and return a float.
 * @return The result of the numerical integration.
 */
double qromo(double (*func)(double), double a, double b, double (*choose)(double(*)(double), double, double, int))
{
	int j;
	double ss,dss,h[JMAXP+1],s[JMAXP+1];

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=(*choose)(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) < EPS*fabs(ss)) return ss;
		}
		s[j+1]=s[j];
		h[j+1]=h[j]/9.0;
	}
	nrerror("Too many steps in routing qromo");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K
