#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "emissionfunctions.h"

double uniform(double a, double b)	{
	double range = b-a;
	return a + (rand()*range)/RAND_MAX;
}

struct vars transform(double t, double phi_k, double k)	{
	struct vars transedVars;

	transedVars.u2 = k/(1-k);
	transedVars.u = pow(transedVars.u2, 0.5);
	transedVars.phi = phi_k*2.*M_PI;
	transedVars.z = t;
	transedVars.z1 = 1-pow(t,2); 
	transedVars.z2 = pow(t,2);

	return transedVars;
}

double ua2(double z, double phi, double u, double u2)	{
	return 1+2.*pow((1-z)/z,0.5)*u*cos(phi)+(1-z)/z*u2; 
}	

double ub2(double z, double phi, double u, double u2)	{
	return 1-2.*pow(z/(1-z),0.5)*u*cos(phi)+z/(1-z)*u2;  
}	

double fc_x(double z, double phi, double u, double u2, double x)	{
	return z*pow(ua2(z,phi,u,u2),1-x/2.)+(1-z)*pow(ub2(z,phi,u,u2),1-x/2.);
}

double fcCorrection(double z, double phi, double u, double u2, double x)	{
	return log(fc_x(z,phi,u,u2,x)); 
}

double nfEquation(double x, double k, double u, double u2, double phi, double z, double correl)	{

	double Hq = 1-((z*(1-z))/(1+u2))*pow(2.*cos(phi)+((1-2.*z)*u)/pow(z*(1-z),0.5),2);

	return (1./(2.*k))*Hq*correl;
}

double caEquation(double x, double k, double t, double u, double u2, double phi, double z, double z1, double z2, double correl, double correlZ1, double correlZ2)	{ 

	double Hg = -4+z*(1-z)/(1+u2)*pow(2.*cos(phi)+((1-2.*z)*u)/pow(z*(1-z),0.5),2);  
	double zCompOne = (1./(1-z1))*(1./2.+1./2.*(1-(1-z1)*u2/z1)/ua2(z1,phi,u,u2)+(1-z1*u2/(1-z1))/ub2(z1,phi,u,u2));
	double zCompTwo = (1./z2)*(1./2.+1./2.*(1-z2*u2/(1-z2))/ub2(z2,phi,u,u2)+(1-(1-z2)*u2/z2)/ua2(z2,phi,u,u2));

	return (1./(2.*k))*(Hg*correl+zCompOne*2.*t*correlZ1+zCompTwo*2.*t*correlZ2);
}

double equation(const char *flag, double x, double k, double phi_k, double t, double u, double u2, double phi, double z, double z1, double z2, double correl, double correlZ1, double correlZ2)	{
	if (strcmp(flag, "nf") == 0)
		return nfEquation(x, k, u, u2, phi, z, correl);
	if (strcmp(flag, "ca") == 0)
		return caEquation(x, k, t, u, u2, phi, z, z1, z2, correl, correlZ1, correlZ2);
	else	{
		printf("invalid flag\n");
		return 0;
	}
}

void fcorrelMC(const char *flag, double x, unsigned long long N, double *I, double *stddev)	{

	*I = 0;
	double r, I2 = 0;
	unsigned long long i;
	double k, t, phi_k, fc, fcZ1, fcZ2;
	double u, u2, phi, z, z1, z2;
	struct vars transedVars;
	srand(time(0)); //seed RNG

	for (i=0; i < N; ++i)	{
		k = uniform(0,1);
		t = uniform(0,1);
		phi_k = uniform(0,1);
		transedVars = transform(t,phi_k,k);
		u = transedVars.u;
		u2 = transedVars.u2;
		phi = transedVars.phi;
		z = transedVars.z;
		z1 = transedVars.z1;
		z2 = transedVars.z2; //TODO integrate ability to reject random number candidates
		fc = fcCorrection(u,u2,phi,z,x);
		fcZ1 = fcCorrection(u,u2,phi,z1,x);
		fcZ2 = fcCorrection(u,u2,phi,z2,x);
		while (isnan(r = equation(flag, x, k, phi_k, t, u, u2, phi, z, z1, z2, fc, fcZ1, fcZ2)))	{
			k = uniform(0,1);
			t = uniform(0,1);
			phi_k = uniform(0,1);
			transedVars = transform(t,phi_k,k);
			u = transedVars.u;
			u2 = transedVars.u2;
			phi = transedVars.phi;
			z = transedVars.z;
			z1 = transedVars.z1;
			z2 = transedVars.z2;
			fc = fcCorrection(u,u2,phi,z,x);
			fcZ1 = fcCorrection(u,u2,phi,z1,x);
			fcZ2 = fcCorrection(u,u2,phi,z2,x);
		}
		*I+=r/N;
		I2+=pow(r,2)/N;
	}

	double mu = *I/N;
	*stddev = pow((1./N)*(I2-2*mu*(*I)+N*pow(mu,2))/N, 0.5);
}

