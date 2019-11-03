#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

double uniform(double a, double b)	{
	double range = b-a;
	return a + (rand()*range)/RAND_MAX;
}

double nfEquation(double x, double k, double z, double phi_k)	{

	while (k == 1)
		k = uniform(0, 1);
	double u = pow(k/(1-k), 0.5);
	double u2 = pow(u,2);
	double phi = phi_k*2*M_PI;

	double fc = z*pow(1+2*pow((1-z)/z,0.5)*u*cos(phi)+((1-z)/z)*u2, 1-x/2)+(1-z)*pow(1-2*pow(z/(1-z),0.5)*u*cos(phi)+(z/(1-z))*u2, 1-x/2); 
	double Hq = 1-((z*(1-z))/(1+u2))*pow(2*cos(phi)+(1-2*z)*u/pow(z*(1-z),0.5),2);

	return (1/(u2*(1+u2)))*(1/(4*M_PI))*Hq*log(fc);
}

double monteCarlo(double (*equation)(double, double, double, double), double x, unsigned long long N)	{
	double I = 0;
	double V = pow(M_PI, 2);
	unsigned long long i;
	double phi_k, k, z;	

	srand(time(0)); //seed RNG

	for (i=0; i < N; ++i)	{
		phi_k = uniform(0, 1);
		z = uniform(0, 1);
		k = uniform(0, 1);
		I += equation(x, k, z, phi_k);
	}

	return (V*I)/N;
}

/* generate all variables between 0 and 1 then transform for substituion equation
 * make sure to check against x=0 analyic result (pi^2/18)
 * remember epislon = 0
 * split up 2S in partial fractions then add coeffs with Hg
 */

int main()	{
	return 0;
}
