// FORSYTHE.h
// включаемый файл, для файла zerion.cpp



#include <math.h>
#include <stdlib.h>
#define REAL double
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define EXP(x) exp(x)
#define SQRT(x) sqrt(x)
#define SIN(x) sin(x)
#define COS(x) cos(x)
#define TAN(x) tan(x)
#define ATAN(x) atan(x)

//enum bool{false, true}; //remove if using MS VC


inline double absval(double x)
{
	return fabs(x);
}
inline int absval(int x)
{
	return abs(x);
}

#include <float.h>

#define FT_FLOAT        0       //'float'
#define FT_DOUBLE       1       //'double'
#define FT_LDOUBLE      2       //'long double'

//Здесь можно переопределить вещественный тип, который будет использоваться в вычислениях.
//Если важна скорость - используйте 'float' (т.е. FT_FLOAT), если нужна точность - используйте 'double' (FT_DOUBLE)
//или 'long double' (FT_LDOUBLE).
#define FLOATTYPE       FT_DOUBLE

#if FLOATTYPE == FT_FLOAT
typedef float Float;
#define EPSILON FLT_EPSILON
#define MAXIMUM FLT_MAX
#elif FLOATTYPE == FT_DOUBLE
typedef double Float;
#define EPSILON DBL_EPSILON
#define MAXIMUM DBL_MAX
#elif FLOATTYPE == FT_LDOUBLE
typedef long double Float;
#define EPSILON LDBL_EPSILON
#define MAXIMUM LDBL_MAX
#else
#error Invalid Float Type
#endif

#define SIGN3(a,b) (b<0) ? (-absval(a)) : absval(a)
#define ABS(x) (((x) >= 0) ? (x) : -(x))
#define SIGN0(x) (((x) > 0) ? 1 : -1)
#define SIGN(x) (((x) == 0) ? 0 : SIGN0(x))
#define SIGN2(a, b) (SIGN(b)*ABS(a))
