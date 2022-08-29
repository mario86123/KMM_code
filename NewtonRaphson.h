#ifndef _newton
#define _newton

#define MAX_IT 100
#define ACC 0.001

/*
 *
 */
double newton(double x0, double x2, double (* f)(double, void *), double (* fdev)(double, void *), void * p);

/*
 * Estimates parameter x of function f in the range [x1,x2]  by means of running the newton iterative algorithm.
 */
double Newton(double x1, double x2, double (* f)(double, void *), double (* fdev)(double, void *), void * p);

#endif
