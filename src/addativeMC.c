#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "emissionfunctions.h"
#include "tests.h"

double fc(double u, double u2, double phi, double z, double x)	{
	return pow(1+u2,(x-1)/2)*(z*pow(ua2(z,phi,u,u2),1-x/2.)+(1-z)*pow(ub2(z,phi,u,u2),-x/2.));
}

double FcorrelMin(double fcorrel, double Cab, double zetaV, double zetaSum, double R)	{
	return 1./zetaV*Cab*(exp(-R*log(zetaV*fcorrel+1+zetaSum))-exp(-R*log(zetaV+1+zetaSum)));
}

double FcorrelMaj(double fcorrel, double Cab, double zetaSum, double R)	{
	return Cab*(exp(-R*log(fcorrel+zetaSum))-exp(-R*log(1+zetaSum)));
}

double fcVsc(double zetaV, double zetaSum, double x, double R, double u, double u2, double phi, double z)	{
	return 1/zetaV*(exp(-R*log(zetaV*fc(u,u2,phi,z,x)+1+zetaSum))-exp(-R*log(zetaV+1+zetaSum)))+(exp(-R*log(fc(u,u2,phi,z,x)+zetaSum))-exp(-R*log(zetaV+1+zetaSum))); 
}

void totalMC(double x, double epsi, double R, unsigned long long N, double *I, double *stddev)	{

	unsigned long long i;
	double phi_k, k, t, zetaV, fcorrel, zetaSum, Cab, r, I2=0;
	double u,u2,phi,z,z1,z2;
	double *transedVars;
	*I = 0;
	
	for (i=0; i < N; ++i)	{
		while ((zetaV = uniform(0,1)) == 0 || zetaV == 1);
		while ((phi_k = uniform(0,1)) == 0 || phi_k == 1);
		while ((k = uniform(0,1)) == 0 || k == 1);
		while ((t = uniform(0,1)) == 0 || t == 1);
		transedVars = transform(t, phi_k, k);	
		u = transedVars[0];
		u2 = transedVars[1];
		phi = transedVars[2];
		z = transedVars[3];
		z1 = transedVars[4];
		z2 = transedVars[5];
		zetaSum = iZeta(ZETA_0, epsi, R, 0);
		Cab = caEquation(x,k,t,u,u2,phi,z,z1,z2,fcVsc(zetaV,zetaSum,x,R,u,u2,phi,z),fcVsc(zetaV,zetaSum,x,R,u,u2,phi,z1),fcVsc(zetaV,zetaSum,x,R,u,u2,phi,z2))+nfEquation(x,k,u,u2,phi,z,fcVsc(zetaV,zetaSum,x,R,u,u2,phi,z));
		fcorrel = fc(u,u2,phi,z,x);
		r = FcorrelMin(fcorrel, Cab, zetaV, zetaSum, R) + FcorrelMaj(fcorrel, Cab, zetaSum, R); 
		*I += r/N; 
		I2 += pow(r,2)/N;
	}

	*stddev = pow((1./N)*(I2-pow(*I,2)),0.5);
}



int main()	{
	srand(time(0));
	double I, stddev;
	int i;
	int step = 20;
	double N = 100000;
	
	test(1, 1e-5, 0.2, N, TEST_X, step, caMC, caComp); 

	return 0;
}
