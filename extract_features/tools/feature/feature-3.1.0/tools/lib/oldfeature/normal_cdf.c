/* $Id: normal_cdf.c,v 1.2 2002/05/15 10:19:55 mliang Exp $ */
/* Copyright (c) 2001 Mike Liang.  All rights reserved. */

/* Normal Cumulative Distribution Function */

#include <math.h>

#define SQR(x) ((x)*(x))
#define SLICES 100


/* Simpson's Rule */
double simpson(double(*f)(double), double a, double b, int n) {
    int m = n/2;
    double h = (b-a)/(2*m);
    double Sn;
    int k;

    Sn = f(a) + 4*f(b-h) + f(b);
    for(k = 1; k < m; k++) {
        Sn += 2*f(a+2*k*h);
        Sn += 4*f(a+2*k*h-h);
    }
    Sn *= h/3.0;

    return Sn;
}


/* Normal Distribution */
inline double normal(double x,double u,double v) {
    return exp(-SQR(x-u)/(2.0*v))/sqrt(v*2.0*M_PI);
}


/* N(0,1) */
inline double normal_pdf(double x) {
    return exp(-SQR(x)/(2.0))/sqrt(2.0*M_PI);
}


/* Cumulative Normal Distribution (1-sided) */
double normal_cdf(double z) {
    if (z < 0)
        return 0.5 - simpson(normal_pdf,0,-z,SLICES);
    return 0.5 + simpson(normal_pdf,0,z,SLICES);
}
