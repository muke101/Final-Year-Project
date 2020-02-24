#ifndef MONTECARLOS_H_
#define MONTECARLOS_H_

double fcorrelMC(double (*equation)(double, double, double, double),     double x, unsigned long long N);

double nfEquation(double x, double k, double t, double phi_k);

double caEquation(double x, double k, double t, double phi_k);

double multiGluonMC(double epsi, double R, unsigned long long N);

#endif
