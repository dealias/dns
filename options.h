#ifndef __options_h__
#define __options_h__ 1

#define IMPLICIT 0
#define COMPLEX 0
#define NUCOMPLEX 0
#define MCREAL 1
#define DOUBLE_PRECISION 1


#include "utils.h"

const Complex I(0.0,1.0);

#if(COMPLEX)
typedef Complex Var;
inline Complex rand_gauss() {return crand_gauss();}
#else
typedef Real Var;
inline Real rand_gauss() {return drand_gauss();}
#endif

#if(NUCOMPLEX)
typedef Var Nu;
#else
typedef Real Nu;
#endif

#endif
