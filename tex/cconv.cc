using namespace std;
#include "Complex.h"
#include "convolution.h"

#define TEST yes
// g++ -g -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse cconv.cc fftw++.cc -lfftw3 -march=native

// icpc -O3 -ansi-alias -malign-double -fp-model fast=2 cconv.cc fftw++.cc -lfftw3

// FFTW: 
// configure --enable-sse2 CC=icpc CFLAGS="-O3 -ansi-alias -malign-double -fp-model fast=2"

// usage: aout [int m] [int pad]

// Number of iterations.
#ifdef TEST
unsigned int N=1;
static const double E=exp(1.0);
const  Complex I(0,1);
#else
unsigned int N=10000000;
#endif
  

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

#ifdef TEST
  for(unsigned int i=0; i < m; i++) f[i]=g[i]=pow(E,i*I);
#else
  for(unsigned int i=0; i < m; i++) f[i]=Complex(3.0,2.0);
  for(unsigned int i=0; i < m; i++) g[i]=Complex(5.0,3.0);
#endif 

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
  for(log2n=0; n > ((unsigned int) 1 << log2n); log2n++);
  n=1 << log2n;
  cout << "n=" << n << endl;
  cout << "m=" << m << endl;
  
  N=N/n;
  if(N < 10) N=10;
  cout << "N=" << N << endl;
  
  int np=pad ? n : m;
  Complex *f=FFTWComplex(np);
  Complex *g=FFTWComplex(np);
#ifdef TEST  
  Complex pseudoh[m];
#endif

  double offset=0.0;
  seconds();
  for(unsigned int i=0; i < N; ++i) {
    seconds();
    offset += seconds();
  }

  if(false) 
    {
    unsigned int L=m;
    Complex *f=FFTWComplex(2*L);
    for(unsigned int i=0; i < L; i++) f[2*i]=d[i];
    for(unsigned int i=0; i < L; i++) f[2*i+1]=-d[i];
    int m=(L+1)/2;
    Complex *u=FFTWComplex(2*(m+1));
    fft0pad fft(m,2,2,f);
    fft.backwards(f,u);
    fft.forwards(f,u);
    FFTWdelete(u);

    for(int i=0; i < 2*m-1; ++i)
      cout << f[2*i] << endl;
    cout << endl;
    for(int i=0; i < 2*m-1; ++i)
      cout << -f[2*i+1] << endl;
  }
  
  double sum=0.0;
  if(!pad) {
    Complex *u=FFTWComplex(m);
    Complex *v=FFTWComplex(m);
    ImplicitConvolution C(m,u,v);
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
    ExplicitConvolution C(n,m,f);
    for(unsigned int i=0; i < N; ++i) {
      // FFTW out-of-place cr routines destroy the input arrays.
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
    cout << endl;
#ifdef TEST    
    for(unsigned int i=0; i < m; i++) pseudoh[i]=f[i];
#endif
  }
  
  if(false)
    {
    DirectConvolution C(m);
    init(f,g);
    Complex *h=FFTWComplex(n);
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
    FFTWdelete(h);
    }
    
#ifdef TEST
  Complex *h=FFTWComplex(n);
  // test accuracy of convolution methods:
  double error=0.0;
  cout << endl;
  for(unsigned int i=0; i < m; i++) {
    // exact solution for test case.
    h[i]=(i+1)*pow(E,i*I);
    cout << h[i] << endl;
  }
  for(unsigned int i=0; i < m; i++) 
    error += abs2(h[i]-pseudoh[i]);
  error=sqrt(error/m);
  cout << "error="<<error<<endl;
  if (error > 1e-8)
    cerr << "Caution! error="<<error<<endl;
  FFTWdelete(h);
#endif    

  
  FFTWdelete(g);
  FFTWdelete(f);
}

