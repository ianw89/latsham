#define FUNC(x) ((*func)(x))

/**
 * @brief Computes the midpoint approximation of the integral of a function.
 *
 * This function calculates the midpoint approximation of the integral of a 
 * given function `func` over the interval [a, b] using `n` subintervals.
 *
 * @param func A pointer to the function to be integrated. The function should 
 *             take a single double argument and return a double.
 * @param a The lower bound of the interval.
 * @param b The upper bound of the interval.
 * @param n The number of subintervals to use for the approximation.
 * @return The midpoint approximation of the integral of the function over the 
 *         interval [a, b].
 */
double midpnt(double (*func)(double), double a, double b, int n)
{
	double x,tnm,sum,del,ddel;
	static double s;
	int it,j;

	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}
#undef FUNC
