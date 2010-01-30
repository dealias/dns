using namespace std;
#include "Complex.h"
#include "convolution.h"

// For FFTW_NO_SIMD:
// g++ -g -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math cconv.cc fftw++.cc -lfftw3 -march=native
//
// Without FFTW_NO_SIMD:
// g++ -g -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse cconv.cc fftw++.cc -lfftw3 -march=native
//
//
// usage: aout [int m]
// optionally specifies the size of m.

// Number of iterations.
unsigned int N=1000;
  
Complex d[]={Complex(-5,3),Complex(3,1),Complex(4,-2),Complex(-3,1),Complex(0,-2),Complex(0,1),Complex(4,0),Complex(-3,-1),Complex(1,2),Complex(2,1),Complex(3,1)};

unsigned int m=sizeof(d)/sizeof(Complex);
unsigned int n=2*m;

using namespace std;

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

inline void init(Complex *f, Complex *g) 
{
//  for(unsigned int i=0; i < m; i++) f[i]=d[i];
//  for(unsigned int i=0; i < m; i++) g[i]=d[i];
  for(unsigned int i=0; i < m; i++) f[i]=Complex(3.0,2.0);
  for(unsigned int i=0; i < m; i++) g[i]=Complex(5.0,3.0);
}

int main(int argc, char* argv[])
{
  fftw::effort |= FFTW_NO_SIMD;
  
  int pad=0;
  
  if (argc >= 2) {
    m=atoi(argv[1]);
  }
 
  if (argc >= 3) {
    pad=atoi(argv[2]);
  }
  
  n=2*m;
  cout << "min padded buffer=" << n << endl;
  unsigned int log2n;
  // Choose next power of 2 for maximal efficiency.
  for(log2n=0; n > (1 << log2n); log2n++);
  n=1 << log2n;
  cout << "n=" << n << endl;
  cout << "m=" << m << endl;
  
  int np=pad ? n : m;
  Complex *f=FFTWComplex(np);
  Complex *g=FFTWComplex(np);
  Complex *h=f;
#ifdef TEST  
  Complex pseudoh[m];
#endif

  double offset=0.0;
  seconds();
  for(int i=0; i < N; ++i) {
    seconds();
    offset += seconds();
  }

  double sum=0.0;
  if(!pad) {
    Complex *u=FFTWComplex(m);
    Complex *v=FFTWComplex(m);
    cconvolution convolve(m,f);
    for(int i=0; i < N; ++i) {
      init(f,g);
      seconds();
      convolve.unpadded(f,g,u,v);
      sum += seconds();
    }
    
    cout << endl;
    cout << "Unpadded:" << endl;
    cout << (sum-offset)/N << endl;
    cout << endl;
    if(m < 100) 
      for(unsigned int i=0; i < m; i++) cout << h[i] << endl;
    else cout << h[0] << endl;
#ifdef TEST    
    for(unsigned int i=0; i < m; i++) pseudoh[i]=h[i];
#endif    
  }
  
  if(pad) {
    sum=0.0;
    cconvolution Convolve(n,m,f);
    for(int i=0; i < N; ++i) {
      // FFTW out-of-place cr routines destroy the input arrays.
      init(f,g);
      seconds();
      Convolve.fft(h,f,g);
      sum += seconds();
    }
    cout << endl;
    cout << "Padded:" << endl;
    cout << (sum-offset)/N << endl;
    cout << endl;
    if(m < 100) 
      for(unsigned int i=0; i < m; i++) cout << h[i] << endl;
    else cout << h[0] << endl;
    cout << endl;
#ifdef TEST    
    for(unsigned int i=0; i < m; i++) pseudoh[i]=h[i];
#endif
  }
  
  cconvolution convolve(m,f);
  init(f,g);
  h=FFTWComplex(n);
  seconds();
  convolve.direct(h,f,g);
  sum=seconds();
  
  cout << endl;
  cout << "Direct:" << endl;
  cout << sum-offset/N << endl;
  cout << endl;

  if(m < 100) 
    for(unsigned int i=0; i < m; i++) cout << h[i] << endl;
  else cout << h[0] << endl;

  // test accuracy of convolution methods:
  double error=0.0;
#ifdef TEST    
  for(unsigned int i=0; i < m; i++) 
    error += abs2(h[i]-pseudoh[i]);
  cout << "error="<<error<<endl;
  if (error > 1e-12)
    cerr << "Caution! error="<<error<<endl;
#endif    
  
//  FFTWdelete(f);
//  FFTWdelete(g);
}

