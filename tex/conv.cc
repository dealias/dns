using namespace std;
#include "convolution.h"

// Compile with: g++ conv.cc fftw++.cc -lfftw3

using namespace std;

int main()
{	
  unsigned int m=683;
  unsigned int n=(2*m-1)*3;
  if(n % 2 == 1) ++n;
  n /= 2;
  
  unsigned int np=n/2+1;
  Complex *f=FFTWComplex(np);
  Complex *g=FFTWComplex(np);
  Complex *h=FFTWComplex(np);
  Complex *d=FFTWComplex(m);
//  Complex d[]={-5,Complex(3,1),Complex(4,-2),Complex(-3,1),Complex(0,-2),Complex(0,1),Complex(4,0),Complex(-3,-1),Complex(1,2),Complex(2,1),Complex(3,1)};
	
//  unsigned int m=sizeof(d)/sizeof(Complex);

  d[0]=1.0;
  for(unsigned int i=1; i < m; i++) d[i]=Complex(3.0,2.0);
  
  for(unsigned int i=0; i < m; i++) f[i]=d[i];
  for(unsigned int i=0; i < m; i++) g[i]=d[i];
  
  convolution convolve(m);
  for(int i=0; i < 10000; ++i)
    convolve.unpadded(h,f,g);
  cout << h[0] << endl;
  
  convolution Convolve(true,m);
  for(int i=0; i < 10000; ++i) {
    for(unsigned int i=0; i < m; i++) f[i]=d[i];
    for(unsigned int i=0; i < m; i++) g[i]=d[i];
    Convolve.fft(h,f,g);
  }
    
  cout << h[0] << endl;
//  for(unsigned int i=0; i < m; i++) cout << h[i] << endl;
//  cout << endl;
  
  /*
  for(unsigned int i=0; i < m; i++) f[i]=d[i];
  for(unsigned int i=0; i < m; i++) g[i]=d[i];
  convolve.direct(h,f,g);
  for(unsigned int i=0; i < m; i++) cout << h[i] << endl;
  */
  
//  FFTWdelete(f);
//  FFTWdelete(g);
//  FFTWdelete(h);
}

