#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

double uniform(double a, double b)	{
	double range = b-a;
	return a + (rand()*range)/RAND_MAX;
}

double ua2(double z, double phi, double u, double u2)	{
	return 1+2.*pow((1-z)/z,0.5)*u*cos(phi)+(1-z)/z*u2; 
}	

double ub2(double z, double phi, double u, double u2)	{
	return 1-2.*pow(z/(1-z),0.5)*u*cos(phi)+z/(1-z)*u2;  
}	

static double fc(double z, double phi, double u, double u2, double x)	{
	return z*pow(ua2(z,phi,u,u2),1-x/2.)+(1-z)*pow(ub2(z,phi,u,u2),1-x/2.);
}

double nfEquation(double x, double k, double t, double phi_k)	{

	if (k == 1 || k == 0 || t == 0 || t == 1)	return 0;
	double u2 = k/(1-k);
	double u = pow(u2, 0.5);
	double phi = phi_k*2*M_PI;
	double z = t;

	double Hq = 1-((z*(1-z))/(1+u2))*pow(2.*cos(phi)+((1-2.*z)*u)/pow(z*(1-z),0.5),2);

	return (1./(2.*k))*Hq*log(fc(z,phi,u,u2,x));
}

double caEquation(double x, double k, double t, double phi_k)	{

	double u2 = k/(1-k);
	double u = pow(u2, 0.5);
	double phi = phi_k*2.*M_PI;
	double z = t;
	double z1 = 1-pow(t,2);
	double z2 = pow(t,2);

	double Hg = -4+z*(1-z)/(1+u2)*pow(2.*cos(phi)+((1-2.*z)*u)/pow(z*(1-z),0.5),2);  
	double zCompOne = (1./(1-z1))*(1./2.+1./2.*(1-(1-z1)*u2/z1)/ua2(z1,phi,u,u2)+(1-z1*u2/(1-z1))/ub2(z1,phi,u,u2));
	double zCompTwo = (1./z2)*(1./2.+1./2.*(1-z2*u2/(1-z2))/ub2(z2,phi,u,u2)+(1-(1-z2)*u2/z2)/ua2(z2,phi,u,u2));

	return (1./(2.*k))*(Hg*log(fc(z,phi,u,u2,x))+zCompOne*2.*t*log(fc(z1,phi,u,u2,x))+zCompTwo*2.*t*log(fc(z2,phi,u,u2,x)));

}

void fcorrelMC(double (*equation)(double, double, double, double), double x, unsigned long long N, double *I, double *stddev)	{

	*I = 0;
	double r, I2 = 0;
	unsigned long long i;
	double phi_k, k, t;
	srand(time(0)); //seed RNG

	for (i=0; i < N; ++i)	{
		phi_k = uniform(0, 1);
		t = uniform(0, 1);
		k = uniform(0, 1);
		while (isnan(r = equation(x, k, t, phi_k)))	{
			phi_k = uniform(0, 1);
			t = uniform(0, 1);
			k = uniform(0, 1);
		}
		*I+=r;
		I2+=pow(r,2);
	}

	double mu = *I/N;
	*stddev = pow((1./N)*(I2-2*mu*(*I)+N*pow(mu,2))/N, 0.5);
}

