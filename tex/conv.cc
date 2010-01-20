using namespace std;
#include "Complex.h"
#include "convolution.h"

// Compile with:
// g++ -g -O3 -msse2 -march=native -mfpmath=sse conv.cc fftw++.cc -lfftw3
// usage: aout [int m]
// optionally specifies the size of m.

using namespace std;

#include <time.h>
#include <sys/times.h>

void timestamp(bool output=true)
{
  static const double ticktime=1.0/sysconf(_SC_CLK_TCK);
  struct tms buf;

  ::times(&buf);
  static double lasttime=0;
  if(output) cout << (buf.tms_utime-lasttime)*ticktime << endl;
  lasttime=buf.tms_utime;
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
  
  Complex *d=FFTWComplex(m);
  d[0]=1.0;
  for(unsigned int i=1; i < m; i++) d[i]=Complex(3.0,2.0);
  
  if(!pad) {
    convolution convolve(m);
    timestamp(false);
    for(unsigned int i=0; i < m; i++) f[i]=d[i];
    for(unsigned int i=0; i < m; i++) g[i]=d[i];
  
    for(int i=0; i < 1000; ++i) {
      
    for(unsigned int i=0; i < m; i++) f[i]=d[i];
    for(unsigned int i=0; i < m; i++) g[i]=d[i];
      convolve.unpadded(h,f,g);
    }
    
    cout << "Unpadded:" << endl;
    timestamp();
    cout << h[0] << endl;
    if(m < 100) 
      for(unsigned int i=0; i < m; i++) cout << h[i] << endl;
  }

  
  if(pad) {
    cout << endl;
    timestamp(false);
    convolution Convolve(n,m);
    for(int i=0; i < 1000; ++i) {
      // FFTW out-of-place cr routines destroy the input arrays.
      for(unsigned int i=0; i < m; i++) f[i]=d[i];
      for(unsigned int i=0; i < m; i++) g[i]=d[i];
      Convolve.fft(h,f,g);
    }
    cout << "Padded:" << endl;
    timestamp();
    cout << h[0] << endl;
    if(m < 100) 
      for(unsigned int i=0; i < m; i++) cout << h[i] << endl;
    cout << endl;
  }
  
#if 1 
  timestamp(false);
  for(unsigned int i=0; i < m; i++) f[i]=d[i];
  for(unsigned int i=0; i < m; i++) g[i]=d[i];
  convolution convolve(m);
  convolve.direct(h,f,g);
  
  cout << endl;
  cout << "Direct/1000:" << endl;
  timestamp();
  cout << h[0] << endl;
  if(m < 100) 
    for(unsigned int i=0; i < m; i++) cout << h[i] << endl;
#endif  
  
  FFTWdelete(f);
  FFTWdelete(g);
  FFTWdelete(h);
}

