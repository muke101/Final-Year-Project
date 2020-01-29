#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

double *arange(double a, double b, size_t x)	{
	double *arr = malloc(x);	
	double step = (b-a)/x;
	int i;
	
	for (i=0; i < x; i++)	{
		printf("%d\n",i);
		arr[i] = a+step*i;
	}
	arr[i] = b;

	return arr;
}

double uniform(double a, double b)  {
   double range = b-a;
   return a + (rand()*range)/RAND_MAX;
}

unsigned fact(unsigned n, unsigned t)	{
	if (n == 1)	{
		return t;	
	}
	return fact(n-1,t*n);	
}

unsigned n = 0; //defining total iZeta calls globally so function can be tail call optimized

void iZeta(double zetaPrev, double epsi, double R, double *res)	{
	double zeta = zetaPrev*pow(uniform(epsi,1),1/R);	
	if (zeta > epsi)	{
		res[0] *= 1./zeta;
		res[1] += zeta;
		n+=1;
		iZeta(zeta, epsi, R, res);	
	}
}

double monteCarlo(double epsi, double R, unsigned long long N)	{
	double I = 0, I2 = 0;
	double res[2];
	unsigned long long i;
	srand(time(0)); //seed RNG

	res[0] = 1; //holds product sum of 1/zeta_i
	for (i=0; i < N; ++i)	{
		double currentProductSum = res[0];
		res[1] = 0; //holds sum of zeta_i
		iZeta(1, epsi, R, res);
		if (res[1] < 1)	{
			I+=res[1];
			//I2+=
		}
		else	{
			res[0] = currentProductSum; //undo invalidated changes to product sum
		}
	}

	I*=res[0]*pow(epsi,R)*pow(R,n)/(N*fact(n,1));
	return I;
}

int main(int argc, char **argv)	{
	int i;
	size_t x = 20;
	unsigned long long N = 10000;
	double epsi = 10e-3;
	double *R = arange(0.1,2,x);
	
	for(i=0; i < x+1; i++)	{
		printf("%f\n",monteCarlo(epsi, R[i], N));	
	}

	return 0;
}
