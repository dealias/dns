using namespace std;
#include "Complex.h"
#include "convolution.h"
#include "utils.h"

// g++ -g -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse conv.cc fftw++.cc -lfftw3 -march=native

// icpc -O3 -ansi-alias -malign-double -fp-model fast=2 conv.cc fftw++.cc -lfftw3

// FFTW: 
// configure --enable-sse2 CC=icpc CFLAGS="-O3 -ansi-alias -malign-double -fp-model fast=2"

// Number of iterations.
unsigned int N0=10000000;
unsigned int N=0;
  
bool Direct=false, Implicit=true, Explicit=false, Test=false;

Complex d[]={-5,Complex(3,1),Complex(4,-2),Complex(-3,1),Complex(0,-2),Complex(0,1),Complex(4,1),Complex(-3,-1),Complex(1,2),Complex(2,1),Complex(3,1),3};

unsigned int m=sizeof(d)/sizeof(Complex);
	
inline void init(Complex *e, Complex *f, Complex *g) 
{
//  for(unsigned int i=0; i < m; i++) f[i]=d[i];
//  for(unsigned int i=0; i < m; i++) g[i]=d[i];
  e[0]=-2.0;
  for(unsigned int i=1; i < m; i++) e[i]=Complex(2.0,-1.0);
  f[0]=1.0;
  for(unsigned int i=1; i < m; i++) f[i]=Complex(3.0,2.0);
  g[0]=2.0;
  for(unsigned int i=1; i < m; i++) g[i]=Complex(5.0,3.0);
}

int main(int argc, char* argv[])
{
#ifndef __SSE2__
  fftw::effort |= FFTW_NO_SIMD;
#endif  
  
#ifdef __GNUC__	
  optind=0;
#endif	
  for (;;) {
    int c = getopt(argc,argv,"deiptN:m:n:");
    if (c == -1) break;
		
    switch (c) {
      case 0:
        break;
      case 'd':
        Direct=true;
        break;
      case 'e':
        Explicit=true;
        Implicit=false;
        break;
      case 'i':
        Implicit=true;
        Explicit=false;
        break;
      case 'p':
        break;
      case 'N':
        N=atoi(optarg);
        break;
      case 't':
        Test=true;
        break;
      case 'm':
        m=atoi(optarg);
        break;
      case 'n':
        N0=atoi(optarg);
        break;
    }
  }

  unsigned int n=4*m-3;
  cout << "min padded buffer=" << n << endl;
  unsigned int log2n;
  // Choose next power of 2 for maximal efficiency.
  for(log2n=0; n > ((unsigned int) 1 << log2n); log2n++);
  n=1 << log2n;
  cout << "n=" << n << endl;
  cout << "m=" << m << endl;
  
  if(N == 0) {
    N=N0/n;
    if(N < 10) N=10;
  }
  cout << "N=" << N << endl;
  
  unsigned int np=Explicit ? n/2+1 : m+1;
    
  Complex *e=ComplexAlign(np);
  Complex *f=ComplexAlign(np);
  Complex *g=ComplexAlign(np);

  Complex *h0=NULL;
  if(Test) h0=ComplexAlign(m);

  
  double offset=0.0, mean=0.0, sigma=0.0;
  double *T=new double[N];
  offset=emptytime(T,N);

  if(Implicit) {
    ImplicitHBiConvolution C(m);
    for(unsigned int i=0; i < N; ++i) {
      init(e,f,g);
      seconds();
      C.convolve(e,f,g);
      T[i]=seconds();
    }
    
    timings(T,N,offset,mean,sigma);
    cout << "\nImplicit:\n" << mean << "\t" << sigma << "\n" << endl;

    if(m < 100) 
      for(unsigned int i=0; i < m; i++) cout << e[i] << endl;
    else cout << e[0] << endl;
    if(Test)
      for(unsigned int i=0; i < m; i++) h0[i]=e[i];
  }
  
  if(Explicit) {
    ExplicitHBiConvolution C(n,m,f);
    for(unsigned int i=0; i < N; ++i) {
      init(e,f,g);
      seconds();
      C.convolve(e,f,g);
      T[i]=seconds();
    }
    
    timings(T,N,offset,mean,sigma);
    cout << "\nExplicit:\n" << mean << "\t" << sigma << "\n" << endl;

    if(m < 100) 
      for(unsigned int i=0; i < m; i++) cout << e[i] << endl;
    else cout << e[0] << endl;
    if(Test) for(unsigned int i=0; i < m; i++) h0[i]=e[i];
  }
  
  if(Direct) {
    DirectHBiConvolution C(m);
    init(e,f,g);
    Complex *h=ComplexAlign(m);
    seconds();
    C.convolve(h,e,f,g);
    mean=seconds();
    
    cout << "\nDirect:\n" << mean << "\t" << "0" << "\n" << endl;

    if(m < 100) 
      for(unsigned int i=0; i < m; i++) cout << h[i] << endl;
    else cout << h[0] << endl;
    deleteAlign(h);
    if(Test) for(unsigned int i=0; i < m; i++) h0[i]=h[i];
  }

  if(Test) {
    Complex *h=ComplexAlign(m);
    double error=0.0;
    cout << endl;
    cout << "Exact:" << endl;
    for(unsigned int i=0; i < m; i++) 
      error += abs2(h[i]-h0[i]);
    cout << "error="<<error<<endl;
    if (error > 1e-12)
      cerr << "Caution! error="<<error<<endl;
    deleteAlign(h);
  }
  
  deleteAlign(g);
  deleteAlign(f);
  deleteAlign(e);
}
