using namespace std;
#include "Complex.h"
#include "convolution.h"

#define TEST 1

// For FFTW_NO_SIMD:
// g++ -g -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math conv.cc fftw++.cc -lfftw3 -march=native
//
// Without FFTW_NO_SIMD:
// g++ -g -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse conv.cc fftw++.cc -lfftw3 -march=native
//
//
// usage: aout [int m]
// optionally specifies the size of m.

// Number of iterations.
unsigned int N=1000;
  
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
  for(unsigned int i=1; i < m-1; i++) f[i]=Complex(3.0,2.0);
  f[m-1]=3.0;
  g[0]=2.0;
  for(unsigned int i=1; i < m-1; i++) g[i]=Complex(5.0,3.0);
  g[m-1]=2.0;
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
  
//  unsigned int m=sizeof(d)/sizeof(Complex);
  unsigned int n=(2*m-1)*3;
  if(n % 2 == 1) ++n;
  n /= 2;
  cout << "min padded buffer=" << n << endl;
  unsigned int log2n;
  // Choose next power of 2 for maximal efficiency.
  for(log2n=0; n > (1 << log2n); log2n++);
  n=1 << log2n;
  cout << "n=" << n << endl;
  cout << "m=" << m << endl;
  
  unsigned int np=pad ? n/2+1 : m;
    
  Complex *f=FFTWComplex(np);
  Complex *g=FFTWComplex(np);
  Complex *h=f;
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
  for(int i=0; i < N; ++i) {
    seconds();
    offset += seconds();
  }

  double sum=0.0;
  if(!pad) {
    unsigned int c=m/2;
    Complex *u=FFTWComplex(c+1);
    Complex *v=FFTWComplex(c+1);
    convolution convolve(m,u);
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
      for(unsigned int i=0; i < m; i++) cout << f[i] << endl;
    else cout << f[0] << endl;
#ifdef TEST    
    for(unsigned int i=0; i < m; i++) pseudoh[i]=f[i];
#endif    
  }
  
  if(pad) {
    sum=0.0;
    convolution Convolve(n,m,f);
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
  
  if(false)
  if(!pad) {
    convolution convolve(m);
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
#ifdef TEST
    double error=0.0;
    cout << endl;
    for(unsigned int i=0; i < m; i++) 
      error += abs2(h[i]-pseudoh[i]);
    cout << "error="<<error<<endl;
    if (error > 1e-12)
      cerr << "Caution! error="<<error<<endl;
#endif
  }
  
//  FFTWdelete(f);
//  FFTWdelete(g);
//  if(!pad) FFTWdelete(h);
}

