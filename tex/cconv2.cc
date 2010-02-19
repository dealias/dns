using namespace std;

#include "Complex.h"
#include "convolution.h"

#include "Array.h"

using namespace Array;

// g++ -g -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse cconv2.cc fftw++.cc -lfftw3 -march=native

// icpc -O3 -ansi-alias -malign-double -fp-model fast=2 cconv2.cc fftw++.cc -lfftw3

// FFTW: 
// configure --enable-sse2 CC=icpc CFLAGS="-O3 -ansi-alias -malign-double -fp-model fast=2"

// usage: aout [int m] [int pad]

// Number of iterations.
unsigned int N=10000000;
unsigned int n=4;
unsigned int m=2;

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

inline void init(array2<Complex>& f, array2<Complex>& g) 
{
  for(unsigned int i=0; i < m; i++) {
    for(unsigned int j=0; j < m; j++) {
      f[i][j]=Complex(3.0,2.0);
      g[i][j]=Complex(5.0,3.0);
//      f[i][j]=i+j;
//      g[i][j]=i+j;
    }
  }
}

int main(int argc, char* argv[])
{
//  fftw::effort |= FFTW_NO_SIMD;
  
  int pad=0;
  
  if(argc >= 2)
    m=atoi(argv[1]);
 
  if(argc >= 3)
    pad=atoi(argv[2]);
  
  n=2*m;
  cout << "min padded buffer=" << n << endl;
  unsigned int log2n;
  // Choose next power of 2 for maximal efficiency.
  for(log2n=0; n > (1 << log2n); log2n++);
  n=1 << log2n;
  cout << "n=" << n << endl;
  cout << "m=" << m << endl;
  
  N=N/(n*n);
  if(N < 10) N=10;
  cout << "N=" << N << endl;
  
  size_t align=sizeof(Complex);
  int np=pad ? n : m;
  array2<Complex> f(np,np,align);
  array2<Complex> g(np,np,align);
#ifdef TEST  
  array2<Complex> pseudoh(m,m,align);
#endif

  double offset=0.0;
  seconds();
  for(int i=0; i < N; ++i) {
    seconds();
    offset += seconds();
  }

  double sum=0.0;
  if(!pad) {
    Complex *u1=FFTWComplex(m);
    Complex *v1=FFTWComplex(m);
    Complex *u2=FFTWComplex(m*m);
    Complex *v2=FFTWComplex(m*m);
    ImplicitConvolution2 C(m,u1,u2);
    for(int i=0; i < N; ++i) {
      init(f,g);
      seconds();
      C.convolve(f,g,u1,v1,u2,v2);
      sum += seconds();
    }
    
    FFTWdelete(u1);
    FFTWdelete(v1);
    FFTWdelete(u2);
    FFTWdelete(v2);
    
    cout << endl;
    cout << "Implicit:" << endl;
    cout << (sum-offset)/N << endl;
    cout << endl;
    if(m < 5) 
      for(unsigned int i=0; i < m; i++) {
        for(unsigned int j=0; j < m; j++)
          cout << f[i][j] << " ";
        cout << endl;
      } else cout << f[0][0] << endl;
    cout << endl;
#ifdef TEST    
    for(unsigned int i=0; i < m; i++) pseudoh[i]=f[i];
#endif    
  }
  
  if(pad) {
    for(int prune=0; prune <= 1; ++prune) {
      sum=0.0;
      ExplicitConvolution2 C(n,m,f,prune);
      for(int i=0; i < N; ++i) {
        init(f,g);
        seconds();
        C.convolve(f,g);
        sum += seconds();
      }
      cout << endl;
      cout << (prune ? "Pruned:" : "Explicit:") << endl;
      cout << (sum-offset)/N << endl;
      cout << endl;
      if(m < 5) 
        for(unsigned int i=0; i < m; i++) {
          for(unsigned int j=0; j < m; j++)
            cout << f[i][j] << "\t";
          cout << endl;
        } else cout << f[0][0] << endl;
    }
#ifdef TEST    
    for(unsigned int i=0; i < m; i++) pseudoh[i]=f[i];
#endif
  }

  if(false)
  {
    array2<Complex> f(m,m,align);
    array2<Complex> g(m,m,align);
    array2<Complex> h(m,m,align);
    DirectConvolution C(m);
    init(f,g);
    seconds();
    C.convolve(h,f,g);
    sum=seconds();
  
    cout << endl;
    cout << "Direct:" << endl;
    cout << sum-offset/N << endl;
    cout << endl;

    if(m < 5) 
      for(unsigned int i=0; i < m; i++) {
        for(unsigned int j=0; j < m; j++)
          cout << h[i][j] << "\t";
        cout << endl;
      } else cout << h[0][0] << endl;

    // test accuracy of convolution methods:
    double error=0.0;
#ifdef TEST    
    for(unsigned int i=0; i < m; i++) 
      error += abs2(h[i]-pseudoh[i]);
    cout << "error="<<error<<endl;
    if (error > 1e-12)
      cerr << "Caution! error="<<error<<endl;
#endif    
  }
}

