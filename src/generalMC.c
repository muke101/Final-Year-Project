#include "emissionfunctions.h"
#define ZETABUFFER 200

double generalZetaI(double *zetaPrev, double epsi, double R, double zetaSum, struct zetaI *leg1, struct zetaI *leg2)	{
	double r;
	*zetaPrev = *zetaPrev*pow(uniform(0,1),1./R);
	if (*zetaPrev > epsi)	{
		r = uniform(0,1);
		if (r < 0.5)	{
			leg1->zeta = *zetaPrev;
			leg1->phi = uniform(0,1)*2.*M_PI;
			return generalZetaI(zetaPrev, epsi, R, zetaSum+(*zetaPrev), ++leg1, leg2);
		}
		else	{
			leg2->zeta = *zetaPrev;
			leg2->phi = uniform(0,1)*2.*M_PI;
			return generalZetaI(zetaPrev, epsi, R, zetaSum+(*zetaPrev), leg1, ++leg2);
		}
	}
	else	{
		return zetaSum;
	}
}

double modIZeta(struct zetaI *leg)	{
	double zeta, phi;
	double sinSum = 0, cosSum = 0;
	int i;

	for (i=0; leg[i].zeta; i++)	{
		zeta = leg[i].zeta;
		phi = leg[i].phi;
		sinSum+=zeta*sin(phi);
	}

	for (i=0; leg[i].zeta; i++)	{
		zeta = leg[i].zeta;
		phi = leg[i].phi;
		cosSum+=zeta*cos(phi);
	}

	return pow(pow(sinSum,2) + pow(cosSum,2),0.5);
}

void generalMC(double epsi, double R, unsigned long long N, double *I, double *stddev)	{
	*I = 0;
	srand(time(0));
	double zetaSum, r, I2 = 0;
	double zeta; 
	unsigned long long i;
	int j;
	struct zetaI leg1[ZETABUFFER];
	struct zetaI leg2[ZETABUFFER];
	for (i=0; i < N; i++)	{
		memset(leg1, 0, ZETABUFFER*sizeof(struct zetaI)); 
		memset(leg2, 0, ZETABUFFER*sizeof(struct zetaI));

		for (j=0; j < ZETABUFFER; j++)	{
			zeta = 1;
			zetaSum = generalZetaI(&zeta, epsi, R, 0, leg1, leg2);	
		}

		r = (zetaSum + modIZeta(leg1) + modIZeta(leg2))/2; 

		*I+=r/N;
		I2+=pow(r,2)/N;
	}

	*stddev = pow((1./N)*(I2-pow(*I,2)),0.5);

}

int main()	{
	double I, stddev;
	generalMC(10e-10, 1, 100000, &I, &stddev);
	printf("I: %.16f, stddev: %.16f\n",I,stddev);
	return 0;
}
