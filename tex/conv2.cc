using namespace std;

#include "Complex.h"
#include "convolution.h"

#include "Array.h"

using namespace Array;

// g++ -g -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse conv2.cc fftw++.cc -lfftw3 -march=native

// icpc -O3 -ansi-alias -malign-double -fp-model fast=2 conv2.cc fftw++.cc -lfftw3

// FFTW: 
// configure --enable-sse2 CC=icpc CFLAGS="-O3 -ansi-alias -malign-double -fp-model fast=2"

// usage: aout [int m] [int pad]

// Number of iterations.
unsigned int N=10000000;
unsigned int nx=0;
unsigned int ny=0;
unsigned int mx=4;
unsigned int my=4;
unsigned int nxp;
unsigned int nyp;

int pad;

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

inline void init(array2<Complex>& f, array2<Complex>& g) 
{
  unsigned int offset=pad ? nx/2-mx+1 : 0;
  unsigned int origin=offset+mx-1;
  unsigned int stop=origin+mx;
  
  for(unsigned int i=offset; i < stop; i++) {
    for(unsigned int j=0; j < my; j++) {
      f[i][j]=Complex(3.0,2.0);
      g[i][j]=Complex(5.0,3.0);
//      f[i][j]=i+j;
//      g[i][j]=i+j;
    }
  }

  f[origin][0]=f[origin][0].re;
  g[origin][0]=g[origin][0].re;

  for(unsigned int i=1; i < mx; i++) {
    f[origin-i][0]=conj(f[origin+i][0]);
    g[origin-i][0]=conj(g[origin+i][0]);
  }
}

unsigned int padding(unsigned int m)
{
  unsigned int n=3*m-2;
  cout << "min padded buffer=" << n << endl;
  unsigned int log2n;
  // Choose next power of 2 for maximal efficiency.
  for(log2n=0; n > ((unsigned int) 1 << log2n); log2n++);
  return 1 << log2n;
}

int main(int argc, char* argv[])
{
//  fftw::effort |= FFTW_NO_SIMD;
  
  pad=0;
  
  if(argc >= 2)
    mx=my=atoi(argv[1]);
 
  if(argc >= 3)
    pad=atoi(argv[2]);
  
  if(argc >= 4)
    my=atoi(argv[3]);
  
  nx=padding(mx);
  ny=padding(my);
  
  cout << "nx=" << nx << ", ny=" << ny << endl;
  cout << "mx=" << mx << ", my=" << my << endl;
  
  N=N/(nx*ny);
  if(N < 10) N=10;
  N=1;
  cout << "N=" << N << endl;
    
  size_t align=sizeof(Complex);
  nxp=pad ? nx : 2*mx-1;
  nyp=pad ? ny/2+1 : my;
  array2<Complex> f(nxp,nyp,align);
  array2<Complex> g(nxp,nyp,align);

  double offset=0.0;
  seconds();
  for(unsigned int i=0; i < N; ++i) {
    seconds();
    offset += seconds();
  }

  double sum=0.0;
  if(!pad) {
    unsigned int c=my/2;
    unsigned int mxpmy=(mx+1)*my;
    Complex *u1=FFTWComplex(c+1);
    Complex *v1=FFTWComplex(c+1);
    Complex *u2=FFTWComplex(mxpmy);
    Complex *v2=FFTWComplex(mxpmy);
    ImplicitHConvolution2 C(mx,my,u1,v1,u2);
    for(unsigned int i=0; i < N; ++i) {
      init(f,g);
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
  
  if(pad) {
    for(unsigned int prune=0; prune <= 1; ++prune) {
      sum=0.0;
      ExplicitHConvolution2 C(nx,ny,mx,my,f,prune);
      for(unsigned int i=0; i < N; ++i) {
        // FFTW out-of-place cr routines destroy the input arrays.
        init(f,g);
        seconds();
        C.convolve(f,g);
        sum += seconds();
      }
      cout << endl;
      cout << (prune ? "Pruned:" : "Explicit:") << endl;
      cout << (sum-offset)/N << endl;
      cout << endl;
      unsigned int offset=nx/2-mx+1;
      if(2*(mx-1)*my < outlimit) 
        for(unsigned int i=offset; i < offset+2*mx-1; i++) {
          for(unsigned int j=0; j < my; j++)
            cout << f[i][j] << "\t";
          cout << endl;
        } else cout << f[offset][0] << endl;
    }
  }
  
//  if(false)
  {
    pad=0;
    unsigned int nxp=2*mx-1;
    array2<Complex> f(nxp,my,align);
    array2<Complex> g(nxp,my,align);
    array2<Complex> h(nxp,my,align);
    DirectHConvolution2 C(mx,my);
    init(f,g);
    seconds();
    C.convolve(h,f,g);
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

