#ifndef MONTECARLOS_H_
#define MONTECARLOS_H_

double uniform(double a, double b);

double ua2(double z, double phi, double u, double u2);

double ub2(double z, double phi, double u, double u2);

double nfEquation(double x, double k, double t, double phi_k);

double caEquation(double x, double k, double t, double phi_k);

double iZeta(double zetaPrev, double epsi, double R, double zetaSum);

double *arange(double a, double b, size_t x);

#endif
