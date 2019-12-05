#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

double uniform(double a, double b)	{
	double range = b-a;
	return a + (rand()*range)/RAND_MAX;
}

double ua2(double z, double phi, double u, double u2)	{
	return 1+2*pow((1-z)/z,0.5)*u*cos(phi)+(1-z)/z*u2; 
}	

double ub2(double z, double phi, double u, double u2)	{
	return 1-2*pow(z/(1-z),0.5)*u*cos(phi)+1/(1-z)*u2;  
}	

double nfEquation(double x, double k, double z, double phi_k)	{

	if (k == 1 || k == 0)	return 0;
	double u2 = k/(1-k);
	double u = pow(u2, 0.5);
	double phi = phi_k*2*M_PI;

	double fc = z*pow(ua2(z,phi,u,u2),1-x/2)+(1-z)*pow(ub2(z,phi,u,u2),1-x/2);
	double Hq = 1-((z*(1-z))/(1+u2))*pow(2*cos(phi)+((1-2*z)*u)/pow(z*(1-z),0.5),2);

	return (1/(2*k))*Hq*log(fc);
}

double caEquation(double x, double k, double z, double phi_k)	{

	if (k == 1)	return 0;
	double u2 = k/(1-k);
	double u = pow(u2, 0.5);
	double phi = phi_k*2*M_PI;
	double z = t;
	double z1 = pow(z,2)+1; 
	double z2 = pow(z,2);

	double fc = z*pow(ua2(z,phi,u,u2),1-x/2)+(1-z)*pow(ub2(z,phi,u,u2),1-x/2);
	double Hg1 = -4+z*(1-z)/(1+u2)*pow(2*cos(phi)+(1-2*z)*u/pow(z*(1-z),0.5),2); 
	double zCompOne = 1/(1-z1)*(1/2+1/2*(1-(1-z1)*u2/z1)/ua2(z1,phi,u,u2)+(1-2*u2/(1-z1)/ub2(z1,phi,u,u2));
	double zCompTwo = 1/z2*(1/2+1/2*(1-2*u2/(1-z))/ub2(z2,phi,u,u2)+(1-(1-z2)*u2/z2)/ua2(z2,phi,u,u2));

	return 2/k*pow(z,0.5),pow(z-1,0.5)*(Hg1+zCompOne+zCompTwo)*log(fc);
}

double stddev;

double monteCarlo(double (*equation)(double, double, double, double), double x, unsigned long long N)	{

	double r, I = 0, I2 = 0;
	unsigned long long i;
	double phi_k, k, z;

	srand(time(0)); //seed RNG

	for (i=0; i < N; ++i)	{
		phi_k = uniform(0, 1);
		z = uniform(0, 1);
		k = uniform(0, 1);
		r = equation(x, k, z, phi_k);
		I+=r;
		I2+=pow(r,2);
	}

	double mu = I/N;
	stddev = pow((I2-2*mu*I+N*pow(mu,2))/N, 0.5);

	return mu; //not *actually* returning mu but works out to be the same computation so saving cycles 
}

double returnStddev()	{
	return stddev;
}

 /* when changing variables, don't need to change the function itself but do need to do the limits and jabobian. You just transform the variables and plug into function. change jacobian to account for d's (measure)
 */

int main()	{
	return 0;
}
