#ifndef EMISSIONFUNCTIONS_H_ 
#define EMISSIONFUNCTIONS_H_

#define ZETA_0 1
#define VARS 6

struct vars	{
	double u;
	double u2;
	double phi;
	double z;
	double z1;
	double z2;
};

double uniform(double a, double b);

double ua2(double z, double phi, double u, double u2);

double ub2(double z, double phi, double u, double u2);

double nfEquation(double x, double k, double u, double u2, double phi, double z, double correl);

double caEquation(double x, double k, double t, double u, double u2, double phi, double z, double z1, double z2, double correl, double correlZ1, double correlZ2);

double iZeta(double zetaPrev, double epsi, double R, double zetaSum);

double fc_x(double z, double u, double u2, double phi, double x);

double fcCorrection(double u, double u2, double phi, double z, double x);

double fcVsc(double zetaV, double zetaSum, double x, double R, double u, double u2, double phi, double z);

struct vars transform(double t, double phi_k, double k);

#endif
