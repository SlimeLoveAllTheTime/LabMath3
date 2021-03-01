#include "Euler.h"
#include <cmath>
#include <cstdlib>

void euler(void (*F)(double t1,double *y1,double *dy1), int n, double h, double *y, double t, double tout)
{
	double *yp = new double[n];
	double *Y = new double[n];
	int i;
	while (t<tout)
	{
		Y[0] = y[0];
		Y[1] = y[1];
		F(t,y,yp);
		for(i=n-1; i>=0; i--) 
			y[i]=Y[i]+h/2*yp[i];
		F(t+h/2,y,yp);
		for(i=n-1; i>=0; i--)
			y[i] = Y[i] + h*yp[i];
		t+=h;
	}
	delete [] yp,Y;
	return;
}