using namespace std;
#include "Complex.h"
#include "convolution.h"

// g++ -g -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse conv.cc fftw++.cc -lfftw3 -march=native

// icpc -O3 -ansi-alias -malign-double -fp-model fast=2 conv.cc fftw++.cc -lfftw3

// FFTW: 
// configure --enable-sse2 CC=icpc CFLAGS="-O3 -ansi-alias -malign-double -fp-model fast=2"

// usage: aout [int m] [int pad]

// Number of iterations.
unsigned int N=10000000;
  
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

Complex d[]={-5,Complex(3,1),Complex(4,-2),Complex(-3,1),Complex(0,-2),Complex(0,1),Complex(4,1),Complex(-3,-1),Complex(1,2),Complex(2,1),Complex(3,1),3};

unsigned int m=sizeof(d)/sizeof(Complex);
	
inline void init(Complex *f, Complex *g) 
{
//  for(unsigned int i=0; i < m; i++) f[i]=d[i];
//  for(unsigned int i=0; i < m; i++) g[i]=d[i];
  f[0]=1.0;
  for(unsigned int i=1; i < m; i++) f[i]=Complex(3.0,2.0);
  g[0]=2.0;
  for(unsigned int i=1; i < m; i++) g[i]=Complex(5.0,3.0);
}

int main(int argc, char* argv[])
{
//  fftw::effort |= FFTW_NO_SIMD;
  
  int pad=0;
  
  if(argc >= 2)
    m=atoi(argv[1]);
 
  if(argc >= 3)
    pad=atoi(argv[2]);
  
  unsigned int n=3*m-2;
  cout << "min padded buffer=" << n << endl;
  unsigned int log2n;
  // Choose next power of 2 for maximal efficiency.
  for(log2n=0; n > ((unsigned int) 1 << log2n); log2n++);
  n=1 << log2n;
  cout << "n=" << n << endl;
  cout << "m=" << m << endl;
  
  N=N/n;
  if(N < 10) N=10;
  cout << "N=" << N << endl;
  
  unsigned int np=pad ? n/2+1 : m;
    
  Complex *f=FFTWComplex(np);
  Complex *g=FFTWComplex(np);
#ifdef TEST  
  Complex pseudoh[m];
#endif

  /*
  Complex *d=FFTWComplex(m);
  d[0]=1.0;
  for(unsigned int i=1; i < m; i++) d[i]=Complex(3.0,2.0);
  */
  
  double offset=0.0;
  seconds();
  for(unsigned int i=0; i < N; ++i) {
    seconds();
    offset += seconds();
  }

  double sum=0.0;
  if(!pad) {
    unsigned int c=m/2;
    Complex *u=FFTWComplex(c+1);
    Complex *v=FFTWComplex(c+1);
    ImplicitHConvolution C(m,u,v);
    for(unsigned int i=0; i < N; ++i) {
      init(f,g);
      seconds();
      C.convolve(f,g,u,v);
      sum += seconds();
    }
    FFTWdelete(u);
    FFTWdelete(v);
    
    cout << endl;
    cout << "Implicit:" << endl;
    cout << (sum-offset)/N << endl;
    cout << endl;
    if(m < 100) 
      for(unsigned int i=0; i < m; i++) cout << f[i] << endl;
    else cout << f[0] << endl;
#ifdef TEST    
    for(unsigned int i=0; i < m; i++) pseudoh[i]=f[i];
#endif    
  }
  
  if(pad) {
    sum=0.0;
    ExplicitHConvolution C(n,m,f);
    for(unsigned int i=0; i < N; ++i) {
      init(f,g);
      seconds();
      C.convolve(f,g);
      sum += seconds();
    }
    
    cout << endl;
    cout << "Explicit:" << endl;
    cout << (sum-offset)/N << endl;
    cout << endl;
    if(m < 100) 
      for(unsigned int i=0; i < m; i++) cout << f[i] << endl;
    else cout << f[0] << endl;
#ifdef TEST    
    for(unsigned int i=0; i < m; i++) pseudoh[i]=f[i];
#endif
  }
  
  if(false)
  {
    DirectHConvolution C(m);
    init(f,g);
    Complex *h=FFTWComplex(m);
    seconds();
    C.convolve(h,f,g);
    sum=seconds();
    
    cout << endl;
    cout << "Direct:" << endl;
    cout << sum-offset/N << endl;
    cout << endl;
    if(m < 100) 
      for(unsigned int i=0; i < m; i++) cout << h[i] << endl;
    else cout << h[0] << endl;

    // test accuracy of convolution methods:
#ifdef TEST
    double error=0.0;
    cout << endl;
    for(unsigned int i=0; i < m; i++) 
      error += abs2(h[i]-pseudoh[i]);
    cout << "error="<<error<<endl;
    if (error > 1e-12)
      cerr << "Caution! error="<<error<<endl;
#endif
    FFTWdelete(h);
  }
  
  FFTWdelete(f);
  FFTWdelete(g);
}
