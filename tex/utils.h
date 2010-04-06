// Utilities for dns/tex

#ifndef __utils_h__
#define __utils_h__ 1
#endif

#include <iostream>
#include <sys/time.h>

inline double seconds()
{
  static timeval lasttime;
  timeval tv;
  gettimeofday(&tv,NULL);
  double seconds=tv.tv_sec-lasttime.tv_sec+
    ((double) tv.tv_usec-lasttime.tv_usec)/1000000.0;
  lasttime=tv;
  return seconds;
}

// timing routines
double emptytime(double *T, unsigned int N)
{
  double val=0.0;
  for(unsigned int i=0; i < N; ++i) {
    seconds();
    T[i]=seconds();
  }
  for(unsigned int i=0; i < N; ++i) 
    val += T[i];
  return val/N;
}

void stdev(double *T, unsigned int N, double mean, double &sigmaL,
           double& sigmaH) 
{
  sigmaL=0.0, sigmaH=0.0;
  unsigned int L=0, H=0;
  for(unsigned int i=0; i < N; ++i) {
    double v=T[i]-mean;
    if(v <= 0) {
      sigmaL += v*v;
      ++L;
    }
    if(v >= 0) {
      ++H;
      sigmaH += v*v;
    }
    
  }
  
  sigmaL=L > 1 ? sqrt(sigmaL/(L-1)) : 0;
  sigmaH=H > 1 ? sqrt(sigmaH/(H-1)) : 0;
}

void timings(const char* text, double *T, unsigned int N)
{
  double sigmaL=0.0, sigmaH=0.0;
  double mean=0.0;
  for(unsigned int i=0; i < N; ++i)
    mean += T[i];
  mean /= N;
  stdev(T,N,mean,sigmaL,sigmaH);
  mean -= emptytime(T,N);
  std::cout << std::endl << text << ":\n" << mean << "\t" << sigmaL << "\t" <<
    sigmaH << std::endl << std::endl;
}
