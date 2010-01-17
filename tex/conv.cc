using namespace std;
#include "Complex.h"
#include "convolution.h"

// Compile with: g++ conv.cc fftw++.cc -lfftw3

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

int main()
{	
//  unsigned int m=8192;
  unsigned int m=5376;
 
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
  
  unsigned int np=n/2+1;
  Complex *f=FFTWComplex(np);
  Complex *g=FFTWComplex(np);
  Complex *h=FFTWComplex(np);
  
  Complex *d=FFTWComplex(m);
  d[0]=1.0;
  for(unsigned int i=1; i < m; i++) d[i]=Complex(3.0,2.0);
  
  for(unsigned int i=0; i < m; i++) f[i]=d[i];
  for(unsigned int i=0; i < m; i++) g[i]=d[i];
  
  convolution convolve(m);
  convolution Convolve(n,m);
  
  timestamp(false);
  
  for(int i=0; i < 1000; ++i)
    convolve.unpadded(h,f,g);
  cout << endl;
  cout << "Unpadded:" << endl;
  timestamp();
  cout << h[0] << endl;
//  for(unsigned int i=0; i < m; i++) cout << h[i] << endl;
  
  for(int i=0; i < 1000; ++i) {
    for(unsigned int i=0; i < m; i++) f[i]=d[i];
    for(unsigned int i=0; i < m; i++) g[i]=d[i];
    Convolve.fft(h,f,g);
  }

  cout << endl;
  cout << "Padded:" << endl;
  timestamp();
  cout << h[0] << endl;
//  for(unsigned int i=0; i < m; i++) cout << h[i] << endl;
  cout << endl;
  
  for(unsigned int i=0; i < m; i++) f[i]=d[i];
  for(unsigned int i=0; i < m; i++) g[i]=d[i];
  convolve.direct(h,f,g);
  
  cout << endl;
  cout << "Direct/1000:" << endl;
  timestamp();
  cout << h[0] << endl;
//  for(unsigned int i=0; i < m; i++) cout << h[i] << endl;
  
  FFTWdelete(f);
  FFTWdelete(g);
  FFTWdelete(h);
}

