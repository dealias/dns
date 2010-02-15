using namespace std;

#include "Complex.h"
#include "convolution.h"

#include "Array.h"

using namespace Array;

// For FFTW_NO_SIMD:
// g++ -g -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math conv2.cc fftw++.cc -lfftw3 -march=native
//
// Without FFTW_NO_SIMD:
// g++ -g -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse conv2.cc fftw++.cc -lfftw3 -march=native
//
//
// usage: aout [int m]
// optionally specifies the size of m.

// Number of iterations.
unsigned int N=100;
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
  
  
//  cout << endl;
//  cout << f << endl;
//  cout << endl;
}

unsigned int padding(unsigned int m)
{
  unsigned int n=(2*m-1)*3;
  if(n % 2 == 1) ++n;
  n /= 2;
  cout << "min padded buffer=" << n << endl;
  unsigned int log2n;
  // Choose next power of 2 for maximal efficiency.
  for(log2n=0; n > (1 << log2n); log2n++);
  return 1 << log2n;
}

int main(int argc, char* argv[])
{
  // Turn off
  fftw::effort |= FFTW_NO_SIMD;
  
  pad=0;
  
  if (argc >= 2) {
    mx=my=atoi(argv[1]);
  }
 
  if (argc >= 3) {
    pad=atoi(argv[2]);
  }
  
  nx=padding(mx);
  ny=padding(my);
  
  cout << "nx=" << nx << ", ny=" << ny << endl;
  cout << "mx=" << mx << ", my=" << my << endl;
  
  size_t align=sizeof(Complex);
  nxp=pad ? nx : 2*mx-1;
  nyp=pad ? ny/2+1 : my;
  array2<Complex> f(nxp,nyp,align);
  array2<Complex> g(nxp,nyp,align);
#ifdef TEST  
  array2<Complex> pseudoh(mx,my,align);
#endif

  double offset=0.0;
  seconds();
  for(int i=0; i < N; ++i) {
    seconds();
    offset += seconds();
  }

  double sum=0.0;
  if(!pad) {
    convolution2 convolve(mx,my,f);
    for(int i=0; i < N; ++i) {
      init(f,g);
      seconds();
      convolve.unpadded(f,g);
      sum += seconds();
    }
    
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
#ifdef TEST    
    for(unsigned int i=0; i < m; i++) pseudoh[i]=f[i];
#endif    
  }
  
  if(pad) {
    for(int prune=0; prune <= 1; ++prune) {
      sum=0.0;
      convolution2 Convolve(nx,ny,mx,my,f,prune);
      for(int i=0; i < N; ++i) {
        // FFTW out-of-place cr routines destroy the input arrays.
        init(f,g);
        seconds();
        Convolve.fft(f,g);
        sum += seconds();
      }
      cout << endl;
      cout << (prune ? "Pruned:" : "Explicit:") << endl;
      cout << (sum-offset)/N << endl;
      cout << endl;
      unsigned int offset=nx/2-mx+1;
      if(nxp*my < outlimit) 
        for(unsigned int i=offset; i < offset+2*mx-1; i++) {
          for(unsigned int j=0; j < my; j++)
            cout << f[i][j] << "\t";
          cout << endl;
        } else cout << f[offset][0] << endl;
    }
#ifdef TEST    
    for(unsigned int i=0; i < m; i++) pseudoh[i]=f[i];
#endif
  }

  if(false)
  if(!pad)
  {
    array2<Complex> h(nxp,nyp,align);
    convolution2 convolve(nx,ny,mx,my,f);
    init(f,g);
    seconds();
    convolve.direct(h,f,g);
    sum=seconds();
  
    cout << endl;
    cout << "Direct:" << endl;
    cout << sum-offset/N << endl;
    cout << endl;

    if(nxp*nyp < outlimit) 
      for(unsigned int i=0; i < nxp; i++) {
        for(unsigned int j=0; j < my; j++)
          cout << h[i][j] << "\t";
        cout << endl;
      } else cout << h[0][0] << endl;

    // test accuracy of convolution methods:
    double error=0.0;
#ifdef TEST    
    for(unsigned int i=0; i < m; i++) 
      error += abs2(h[i]-pseudoh[i]);
    cout << "error="<<error<<endl;
    if (error > 1e-12)
      cerr << "Caution! error="<<error<<endl;
#endif    
  }
  
  
//  FFTWdelete(f);
//  FFTWdelete(g);
}

