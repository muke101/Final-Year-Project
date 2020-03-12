#include "emissionfunctions.h"
#include "tests.h"

double fc(double u, double u2, double phi, double z, double x)	{
	return pow(1+u2,(x-1)/2)*fc_x(z, phi, u, u2, x);
}

double FcorrelMin(double fcorrel, double zetaV, double zetaSum, double R)	{
	return 1./zetaV*(exp(-R*log(zetaV*fcorrel+1+zetaSum))-exp(-R*log(zetaV+1+zetaSum)));
}

double FcorrelMaj(double fcorrel, double zetaSum, double R)	{
	return 1./R*(exp(-R*log(fcorrel+zetaSum))-exp(-R*log(1+zetaSum)));
}

double fcVsc(double zetaV, double zetaSum, double x, double R, double u, double u2, double phi, double z)	{
	return FcorrelMin(fc(u,u2,phi,z,x),zetaV,zetaSum,R)+FcorrelMaj(fc(u,u2,phi,z,x),zetaSum,R);
}

void totalMC(double x, double epsi, double R, unsigned long long N, double *I, double *stddev)	{

	unsigned long long i;
	double phi_k, k, t, zetaV, fcorrel, zetaSum, Cab, fc, fcZ1, fcZ2, r, I2=0;
	double u,u2,phi,z,z1,z2;
	struct vars transedVars;
	*I = 0;
	
	for (i=0; i < N; ++i)	{
		while ((zetaV = uniform(0,1)) == 0 || zetaV == 1);
		while ((phi_k = uniform(0,1)) == 0 || phi_k == 1);
		while ((k = uniform(0,1)) == 0 || k == 1);
		while ((t = uniform(0,1)) == 0 || t == 1);
		transedVars = transform(t, phi_k, k);	
		u = transedVars.u;
		u2 = transedVars.u2;
		phi = transedVars.phi;
		z = transedVars.z;
		z1 = transedVars.z1;
		z2 = transedVars.z2;
		zetaSum = iZeta(ZETA_0, epsi, R, 0);
		fc = fcVsc(zetaV,zetaSum,x,R,u,u2,phi,z);
		fcZ1 = fcVsc(zetaV,zetaSum,x,R,u,u2,phi,z1);
		fcZ2 = fcVsc(zetaV,zetaSum,x,R,u,u2,phi,z2);
		Cab = caEquation(x,k,t,u,u2,phi,z,z1,z2,fc,fcZ1,fcZ2)+nfEquation(x,k,u,u2,phi,z,fc);
		r = Cab;
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

	test(1, 1e-10, 0.2, N, TEST_X, step, nfMC, nfComp); 

	return 0;
}
