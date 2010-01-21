using namespace std;
#include "Complex.h"
#include "convolution.h"

// Compile with:
// g++ -g -O3 -msse2 -march=native -mfpmath=sse conv.cc fftw++.cc -lfftw3
//
// g++ -g -O3 -fomit-frame-pointer -fstrict-aliasing -ffast-math -march=native -msse2 -mfpmath=sse conv.cc fftw++.cc -lfftw3
//
// g++ -g -O3 -fomit-frame-pointer -fstrict-aliasing -ffast-math -mpentiumpro -msse2 -mfpmath=sse conv.cc fftw++.cc -lfftw3 -I$HOME/include -L$HOME/lib
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

int main(int argc, char* argv[])
{
  unsigned int m=4096;
//  unsigned int m=2731;
  
  int pad=0;
  
  if (argc >= 2) {
    m=atoi(argv[1]);
  }
 
  if (argc >= 3) {
    pad=atoi(argv[2]);
  }
  
//  Complex d[]={-5,Complex(3,1),Complex(4,-2),Complex(-3,1),Complex(0,-2),Complex(0,1),Complex(4,0),Complex(-3,-1),Complex(1,2),Complex(2,1),Complex(3,1)};
	
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
  Complex *h=FFTWComplex(np);
  Complex pseudoh[m];


  Complex *d=FFTWComplex(m);
  d[0]=1.0;
  for(unsigned int i=1; i < m; i++) d[i]=Complex(3.0,2.0);
  
  double offset=0.0;
  seconds();
  for(int i=0; i < N; ++i) {
    seconds();
    offset += seconds();
  }

  double sum=0.0;
  if(!pad) {
    convolution convolve(m);
    for(int i=0; i < N; ++i) {
      for(unsigned int i=0; i < m; i++) f[i]=d[i];
      for(unsigned int i=0; i < m; i++) g[i]=d[i];
      seconds();
      convolve.unpadded(h,f,g);
      sum += seconds();
    }
    
    cout << endl;
    cout << "Unpadded:" << endl;
    cout << (sum-offset)/N << endl;
    if(m < 100) 
      for(unsigned int i=0; i < m; i++) cout << h[i] << endl;
    else cout << h[0] << endl;
    for(unsigned int i=0; i < m; i++) pseudoh[i]=h[i];
  }
  
  if(pad) {
    sum=0.0;
    convolution Convolve(n,m);
    for(int i=0; i < N; ++i) {
      // FFTW out-of-place cr routines destroy the input arrays.
      for(unsigned int i=0; i < m; i++) f[i]=d[i];
      for(unsigned int i=0; i < m; i++) g[i]=d[i];
      seconds();
      Convolve.fft(h,f,g);
      sum += seconds();
    }
    cout << endl;
    cout << "Padded:" << endl;
    cout << (sum-offset)/N << endl;
    if(m < 100) 
      for(unsigned int i=0; i < m; i++) cout << h[i] << endl;
    else cout << h[0] << endl;
    cout << endl;
    for(unsigned int i=0; i < m; i++) pseudoh[i]=h[i];
  }
  
#if 1 
  sum=0.0;
  for(unsigned int i=0; i < m; i++) f[i]=d[i];
  for(unsigned int i=0; i < m; i++) g[i]=d[i];
  convolution convolve(m);
  seconds();
  convolve.direct(h,f,g);
  sum +=seconds();
  
  cout << endl;
  cout << "Direct:" << endl;
  cout << (sum-offset)/N << endl;

  if(m < 100) 
    for(unsigned int i=0; i < m; i++) cout << h[i] << endl;
  else cout << h[0] << endl;

  // test accuracy of convolution methods:
  double error=0.0;
  for(unsigned int i=0; i < m; i++) 
    error += abs2(h[i]-pseudoh[i]);
  cout << "error="<<error<<endl;
  if (error > 1e-12)
    cerr << "Caution! error="<<error<<endl;
#endif  
  
  FFTWdelete(f);
  FFTWdelete(g);
  FFTWdelete(h);
}

