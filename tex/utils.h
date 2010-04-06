// Utilities for dns/tex

#ifndef __utils_h__
#define __utils_h__ 1
#endif

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
  return val;
}

void rmoffset(double *T, unsigned int N, double offset)
{
  for(unsigned int i=0; i < N; ++i) {
    T[i] -= offset/N;
  }
}

double stdev(double *T, unsigned int N, double mean) 
{
  double sigma=0.0;
  for(unsigned int i=0; i < N; ++i) {
    double v=T[i]-mean;
    sigma += v*v;
  }
  
  return N > 1 ? sqrt(sigma/(N-1)) : 0;
}

void timings(double *T, unsigned int N, double &offset, double &mean, 
	     double &sigma)
{
  mean=0.0;
  rmoffset(T,N,offset);
  for(unsigned int i=0; i < N; ++i) 
    mean += T[i];
  mean /= N;
  sigma=stdev(T,N,mean);
}
