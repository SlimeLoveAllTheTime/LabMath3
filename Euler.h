#ifndef _EULER
#define _EULER

void euler(void(*F)(double t1,double* y1, double* dy1), int n, double h, double * y, double t, double tout);

#endif