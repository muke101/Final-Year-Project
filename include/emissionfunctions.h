#ifndef EMISSIONFUNCTIONS_H_ 
#define EMISSIONFUNCTIONS_H_

#define ZETA_0 1

double uniform(double a, double b);

double ua2(double z, double phi, double u, double u2);

double ub2(double z, double phi, double u, double u2);

double nfEquation(double x, double k, double t, double phi_k);

double caEquation(double x, double k, double t, double phi_k);

double iZeta(double zetaPrev, double epsi, double R, double zetaSum);

double FcorrelMin(double fcorrel, double Cab, double zetaV, double zetaSum, double R);

double FcorrelMaj(double fcorrel, double Cab, double zetaSum, double R);

double fc(double t, double phi_k, double k, double x);

#endif
