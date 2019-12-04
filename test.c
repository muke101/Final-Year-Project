#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

double uniform(double a, double b)	{
	double range = b-a;
	return a + (rand()*range)/RAND_MAX;
}

double nfEquation(double x, double k, double z, double phi_k)	{

	while (k == 1)	k = uniform(0, 1);
	double u2 = k/(1-k);
	double u = pow(u2, 0.5);
	double phi = phi_k*(2*M_PI);

	double fc = z*pow(1+2*pow((1-z)/z,0.5)*u*cos(phi)+((1-z)*u2)/z,(1-x/2))+(1-z)*pow(1-2*pow(z/(1-z),0.5)*u*cos(phi)+(z*u2)/(1-z),(1-x/2));
	double Hq = 1-((z*(1-z))/(1+u2))*pow(2*cos(phi)+(u*(1-2*z))/pow(z*(1-z),0.5),2);

	return 1/(u2*(1+u2))*1/(2*M_PI)*1/2*Hq*log(fc);

}

double monteCarlo(double (*equation)(double, double, double, double), double x, unsigned long long N)	{
	double I = 0;
	unsigned long long i;
	double phi_k, k, z;	

	srand(time(0)); //seed RNG
	
	for (i=0; i < N; i++)	{
		phi_k = uniform(0, 1);
		z = uniform(0, 1);
		k = uniform(0, 1);
		I += equation(x, k, z, phi_k);
	}

	return I/N;
}

int main()	{
	printf("%f\n",monteCarlo(nfEquation,0,10000));
	return 0;
}
