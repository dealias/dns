#include "Complex.h"
#include "convolution.h"

#ifdef __SSE2__
const union uvec sse2_pm = {
  { 0x00000000,0x00000000,0x00000000,0x80000000 }
};

const union uvec sse2_mm = {
  { 0x00000000,0x80000000,0x00000000,0x80000000 }
};
#endif

const double sqrt3=sqrt(3.0);
const double hsqrt3=0.5*sqrt3;

const Complex hSqrt3(hsqrt3,hsqrt3);
const Complex mhsqrt3(-hsqrt3,-hsqrt3);
const Complex mhalf(-0.5,-0.5);
const Complex zeta3(-0.5,hsqrt3);
