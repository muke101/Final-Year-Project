#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

double nfEquation(double x, double k, double z, double phi)	{
	double fc = z*pow(1+2*pow((1-z)/z,0.5)*tan(k)*cos(phi)+pow(tan(k),2)*((1-z)/z),(1-x/2))+(1-z)*pow(1-2*pow(z/(1-z),0.5)*tan(k)*cos(phi)+(z*pow(tan(k),2)/(1-z)), (1-x/2));
	double Hq = 1-((z*(1-z))/(1+pow(tan(k),2)))*pow(2*cos(phi)+tan(k)*(1-2*z)/pow(z*(1-z),0.5), 2);
	return (log(fc)*Hq)/(4*M_PI);
}

double uniform(double a, double b)	{
	double range = b-a;
	return a + (rand()*range)/RAND_MAX;
}

double monteCarlo(double (*equation)(double, double, double, double), double x, unsigned long long N)	{
	double I = 0;
	double V = pow(M_PI, 2);
	unsigned long long i;
	double phi, k, z;	

	srand(time(0)); //seed RNG

	for (i=0; i < N; ++i)	{
		phi = uniform(0, 2*M_PI);
		z = uniform(0, 1);
		while ((k=uniform(0, M_PI/2)) == M_PI/2)
			k = uniform(0, M_PI/2);
		I += equation(x, k, z, phi);
	}

	return (V*I)/N;
}

int main()	{
	return 0;
}
