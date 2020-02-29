#ifndef TESTS_H_
#define TESTS_H_

double *arange(double a, double b, int x);

void caMC(double x, double epsi, double R, unsigned long long N, double *I);

void caComp(double x, double epsi, double R, unsigned long long N, double *I);

void nfMC(double x, double epsi, double R, unsigned long long N, double *I);

void nfComp(double x, double epsi, double R, unsigned long long N, double *I);

void test(double defX, double defEpsi, double defR, unsigned long long N, const char flags, int step, void(*MCfunc)(double, double, double, unsigned long long, double*), void(*CompFunc)(double, double, double, unsigned long long, double*));

#endif
