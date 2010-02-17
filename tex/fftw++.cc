#include "cmult-sse2.h"
#include "fftw++.h"

std::ifstream fftw::ifWisdom;
std::ofstream fftw::ofWisdom;
bool fftw::Wise=false;

// User settings:
unsigned int fftw::effort=FFTW_PATIENT;
const char *fftw::WisdomName="wisdom3.txt";

#ifdef __SSE2__
const union uvec sse2_pm = {
  { 0x00000000,0x00000000,0x00000000,0x80000000 }
};

const union uvec sse2_mm = {
  { 0x00000000,0x80000000,0x00000000,0x80000000 }
};
#endif

double sqrt3=sqrt(3.0);
double hsqrt3=0.5*sqrt3;

Complex hSqrt3(hsqrt3,hsqrt3);
Complex mhsqrt3(-hsqrt3,-hsqrt3);
Complex mhalf(-0.5,-0.5);
Complex zeta3(-0.5,hsqrt3);
Complex one(1.0,0.0);
