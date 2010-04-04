using namespace std;

#include "Complex.h"
#include "convolution.h"

#include "Array.h"

using namespace Array;

// g++ -g -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse conv3.cc fftw++.cc -lfftw3 -march=native

// icpc -O3 -ansi-alias -malign-double -fp-model fast=2 conv3.cc fftw++.cc -lfftw3

// FFTW: 
// configure --enable-sse2 CC=icpc CFLAGS="-O3 -ansi-alias -malign-double -fp-model fast=2"

// Number of iterations.
unsigned int N0=10000000;
unsigned int N=0;
unsigned int nx=0;
unsigned int ny=0;
unsigned int nz=0;
unsigned int mx=4;
unsigned int my=4;
unsigned int mz=4;
unsigned int nxp;
unsigned int nyp;
unsigned int nzp;

bool Direct=false, Implicit=true;

unsigned int outlimit=300;

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

inline void init(array3<Complex>& f, array3<Complex>& g) 
{
  unsigned int xorigin=mx-1;
  unsigned int yorigin=my-1;
  unsigned int xstop=xorigin+mx;
  unsigned int ystop=yorigin+my;
  
  for(unsigned int i=0; i < xstop; ++i) {
    for(unsigned int j=0; j < ystop; ++j) {
      for(unsigned int k=0; k < mz; ++k) {
        f[i][j][k]=Complex(3.0,2.0);
        g[i][j][k]=Complex(5.0,3.0);
      }
    }
  }

  f[xorigin][yorigin][0]=f[xorigin][yorigin][0].re;
  g[xorigin][yorigin][0]=g[xorigin][yorigin][0].re;
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
#ifndef __SSE2__
  fftw::effort |= FFTW_NO_SIMD;
#endif  
  
#ifdef __GNUC__	
  optind=0;
#endif	
  for (;;) {
    int c = getopt(argc,argv,"deiptN:m:x:y:");
    if (c == -1) break;
		
    switch (c) {
      case 0:
        break;
      case 'd':
        Direct=true;
        break;
      case 'e':
        Implicit=false;
        break;
      case 'i':
        Implicit=true;
        break;
      case 'p':
        Implicit=false;
        break;
      case 'N':
        N=atoi(optarg);
        break;
      case 'm':
        mx=my=mz=atoi(optarg);
        break;
      case 'x':
        mx=atoi(optarg);
        break;
      case 'y':
        my=atoi(optarg);
        break;
      case 'z':
        mz=atoi(optarg);
        break;
    }
  }

  nx=padding(mx);
  ny=padding(my);
  nz=padding(mz);
  
  cout << "nx=" << nx << ", ny=" << ny << ", nz=" << ny << endl;
  cout << "mx=" << mx << ", my=" << my << ", mz=" << mz << endl;
  
  if(N == 0) {
    N=N0/(nx*ny*nz);
    if(N < 10) N=10;
  }
  cout << "N=" << N << endl;
    
  size_t align=sizeof(Complex);
  nxp=2*mx-1;
  nyp=2*my-1;
  nzp=mz;
  array3<Complex> f(nxp,nyp,nzp,align);
  array3<Complex> g(nxp,nyp,nzp,align);

  double offset=0.0;
  seconds();
  for(unsigned int i=0; i < N; ++i) {
    seconds();
    offset += seconds();
  }

  double sum=0.0;
  if(Implicit) {
    ImplicitHConvolution3 C(mx,my,mz);
    for(unsigned int i=0; i < N; ++i) {
      init(f,g);
      seconds();
      C.convolve(f,g);
      sum += seconds();
    }
    
    cout << endl;
    cout << "Implicit:" << endl;
    cout << (sum-offset)/N << endl;
    cout << endl;
    if(nxp*nyp*mz < outlimit) {
      for(unsigned int i=0; i < nxp; ++i) {
        for(unsigned int j=0; j < nyp; ++j) {
          for(unsigned int k=0; k < mz; ++k)
            cout << f[i][j][k] << "\t";
          cout << endl;
        }
        cout << endl;
      }
    } else cout << f[0][0][0] << endl;
  }
  
  if(Direct) {
    unsigned int nxp=2*mx-1;
    unsigned int nyp=2*my-1;
    array3<Complex> h(nxp,nyp,mz,align);
    array3<Complex> f(nxp,nyp,mz,align);
    array3<Complex> g(nxp,nyp,mz,align);
    DirectHConvolution3 C(mx,my,mz);
    init(f,g);
    seconds();
    C.convolve(h,f,g);
    sum=seconds();
  
    cout << endl;
    cout << "Direct:" << endl;
    cout << sum-offset/N << endl;
    cout << endl;

    if(nxp*nyp*mz < outlimit)
      for(unsigned int i=0; i < nxp; ++i) {
        for(unsigned int j=0; j < nyp; ++j) {
          for(unsigned int k=0; k < mz; ++k)
            cout << h[i][j][k] << "\t";
          cout << endl;
        }
        cout << endl;
      }
    else cout << h[0][0][0] << endl;
    
  }
}
