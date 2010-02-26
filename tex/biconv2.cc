using namespace std;

#include "Complex.h"
#include "convolution.h"

#include "Array.h"

using namespace Array;

// g++ -g -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse conv2.cc fftw++.cc -lfftw3 -march=native

// icpc -O3 -ansi-alias -malign-double -fp-model fast=2 conv2.cc fftw++.cc -lfftw3

// FFTW: 
// configure --enable-sse2 CC=icpc CFLAGS="-O3 -ansi-alias -malign-double -fp-model fast=2"

// Number of iterations.
unsigned int N0=10000000;
unsigned int N=0;
unsigned int nx=0;
unsigned int ny=0;
unsigned int mx=4;
unsigned int my=4;
unsigned int nxp;
unsigned int nyp;

bool Direct=false, Implicit=true, Explicit=false, Pruned=false;

unsigned int outlimit=100;

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

inline void init(array2<Complex>& e, array2<Complex>& f, array2<Complex>& g) 
{
  unsigned int offset=Explicit ? nx/2-mx+1 : 0;
  unsigned int origin=offset+mx-1;
  unsigned int stop=origin+mx;
  
  for(unsigned int i=offset; i < stop; i++) {
    for(unsigned int j=0; j < my; j++) {
      e[i][j]=Complex(2.0,-1.0);
      f[i][j]=Complex(3.0,2.0);
      g[i][j]=Complex(5.0,3.0);
//      f[i][j]=i+j;
//      g[i][j]=i+j;
    }
  }

  e[origin][0]=e[origin][0].re;
  f[origin][0]=f[origin][0].re;
  g[origin][0]=g[origin][0].re;

  for(unsigned int i=1; i < mx; i++) {
    e[origin-i][0]=conj(e[origin+i][0]);
    f[origin-i][0]=conj(f[origin+i][0]);
    g[origin-i][0]=conj(g[origin+i][0]);
  }
}

unsigned int padding(unsigned int m)
{
  unsigned int n=4*m-3;
  cout << "min padded buffer=" << n << endl;
  unsigned int log2n;
  // Choose next power of 2 for maximal efficiency.
  for(log2n=0; n > ((unsigned int) 1 << log2n); log2n++);
  return 1 << log2n;
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
    int c = getopt(argc,argv,"deiptN:x:y:");
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
        Pruned=false;
        break;
      case 'i':
        Implicit=true;
        Explicit=false;
        break;
      case 'p':
        Explicit=true;
        Implicit=false;
        Pruned=true;
        break;
      case 'N':
        N=atoi(optarg);
        break;
      case 'x':
        mx=atoi(optarg);
        break;
      case 'y':
        my=atoi(optarg);
        break;
    }
  }

  nx=padding(mx);
  ny=padding(my);
  
  cout << "nx=" << nx << ", ny=" << ny << endl;
  cout << "mx=" << mx << ", my=" << my << endl;
  
  if(N == 0) {
    N=N0/(nx*ny);
    if(N < 10) N=10;
  }
  cout << "N=" << N << endl;
    
  size_t align=sizeof(Complex);
  nxp=Explicit ? nx : 2*mx-1;
  nyp=Explicit ? ny/2+1 : my;
  array2<Complex> e(nxp,nyp,align);
  array2<Complex> f(nxp,nyp,align);
  array2<Complex> g(nxp,nyp,align);

  double offset=0.0;
  seconds();
  for(unsigned int i=0; i < N; ++i) {
    seconds();
    offset += seconds();
  }

  double sum=0.0;
  if(Implicit && false) {
    unsigned int c=my/2;
    unsigned int mxpmy=(mx+1)*my;
    Complex *u1=FFTWComplex(c+1);
    Complex *v1=FFTWComplex(c+1);
    Complex *u2=FFTWComplex(mxpmy);
    Complex *v2=FFTWComplex(mxpmy);
    ImplicitHConvolution2 C(mx,my,u1,v1,u2);
    for(unsigned int i=0; i < N; ++i) {
      init(e,f,g);
      seconds();
      C.convolve(f,g,u1,v1,u2,v2);
      sum += seconds();
    }
    
    FFTWdelete(v2);
    FFTWdelete(u2);
    FFTWdelete(v1);
    FFTWdelete(u1);
    
    cout << endl;
    cout << "Implicit:" << endl;
    cout << (sum-offset)/N << endl;
    cout << endl;
    if(nxp*my < outlimit)
      for(unsigned int i=0; i < nxp; i++) {
        for(unsigned int j=0; j < my; j++)
          cout << f[i][j] << " ";
        cout << endl;
      } else cout << f[0][0] << endl;
    cout << endl;
  }
  
  if(Explicit) {
    sum=0.0;
    ExplicitHBiConvolution2 C(nx,ny,mx,my,f,Pruned);
    for(unsigned int i=0; i < N; ++i) {
      init(e,f,g);
      seconds();
      C.convolve(e,f,g);
      sum += seconds();
    }
    cout << endl;
    cout << (Pruned ? "Pruned:" : "Explicit:") << endl;
    cout << (sum-offset)/N << endl;
    cout << endl;
    unsigned int offset=nx/2-mx+1;
    if(2*(mx-1)*my < outlimit) 
      for(unsigned int i=offset; i < offset+2*mx-1; i++) {
        for(unsigned int j=0; j < my; j++)
          cout << e[i][j] << "\t";
        cout << endl;
      } else cout << e[offset][0] << endl;
  }
  
  if(Direct) {
    Explicit=0;
    unsigned int nxp=2*mx-1;
    array2<Complex> e(nxp,my,align);
    array2<Complex> f(nxp,my,align);
    array2<Complex> g(nxp,my,align);
    array2<Complex> h(nxp,my,align);
    DirectHBiConvolution2 C(mx,my);
    init(e,f,g);
    seconds();
    C.convolve(h,e,f,g);
    sum=seconds();
  
    cout << endl;
    cout << "Direct:" << endl;
    cout << sum-offset/N << endl;
    cout << endl;

    if(nxp*my < outlimit)
      for(unsigned int i=0; i < nxp; i++) {
        for(unsigned int j=0; j < my; j++)
          cout << h[i][j] << "\t";
        cout << endl;
      } else cout << h[0][0] << endl;
  }
}