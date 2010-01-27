using namespace std;

#include "Complex.h"
#include "convolution.h"

#include "Array.h"

using namespace Array;

// For FFTW_NO_SIMD:
// g++ -g -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math cconv2.cc fftw++.cc -lfftw3 -march=native
//
// Without FFTW_NO_SIMD:
// g++ -g -O3 -DNDEBUG -fomit-frame-pointer -fstrict-aliasing -ffast-math -msse2 -mfpmath=sse cconv2.cc fftw++.cc -lfftw3 -march=native
//
//
// usage: aout [int m]
// optionally specifies the size of m.

// Number of iterations.
unsigned int N=100;
unsigned int n=4;
unsigned int m=2;

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
  for(unsigned int i=0; i < m; i++) {
    for(unsigned int j=0; j < m; j++) {
      f[i][j]=Complex(3.0,2.0);
      g[i][j]=Complex(5.0,3.0);
//      f[i][j]=i+j;
//      g[i][j]=i+j;
    }
  }
}

int main(int argc, char* argv[])
{
  // Turn off
  fftw::effort |= FFTW_NO_SIMD;
  
  int pad=0;
  
  if (argc >= 2) {
    m=atoi(argv[1]);
  }
 
  if (argc >= 3) {
    pad=atoi(argv[2]);
  }
  
  n=2*m;
  cout << "min padded buffer=" << n << endl;
  unsigned int log2n;
  // Choose next power of 2 for maximal efficiency.
  for(log2n=0; n > (1 << log2n); log2n++);
  n=1 << log2n;
  cout << "n=" << n << endl;
  cout << "m=" << m << endl;
  
  size_t align=sizeof(Complex);
  int np=pad ? n : m;
  array2<Complex> f(np,np,align);
  array2<Complex> g(np,np,align);
  array2<Complex> h(np,np,align);
#ifdef TEST  
  array2<Complex> pseudoh(m,m,align);
#endif

  double offset=0.0;
  seconds();
  for(int i=0; i < N; ++i) {
    seconds();
    offset += seconds();
  }

  double sum=0.0;
  if(!pad) {
    cconvolution2 convolve(m,f);
    for(int i=0; i < N; ++i) {
      init(f,g);
      seconds();
      convolve.unpadded(f,g);
      sum += seconds();
    }
    
    cout << endl;
    cout << "Unpadded:" << endl;
    cout << (sum-offset)/N << endl;
    cout << endl;
    if(m < 5) 
      for(unsigned int i=0; i < m; i++) {
        for(unsigned int j=0; j < m; j++)
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
      cconvolution2 Convolve(n,m,f,prune);
      for(int i=0; i < N; ++i) {
        // FFTW out-of-place cr routines destroy the input arrays.
        init(f,g);
        seconds();
        Convolve.fft(f,g);
        sum += seconds();
      }
      cout << endl;
      cout << (prune ? "Pruned:" : "Padded:") << endl;
      cout << (sum-offset)/N << endl;
      cout << endl;
      if(m < 5) 
        for(unsigned int i=0; i < m; i++) {
          for(unsigned int j=0; j < m; j++)
            cout << f[i][j] << "\t";
          cout << endl;
        } else cout << f[0][0] << endl;
    }
#ifdef TEST    
    for(unsigned int i=0; i < m; i++) pseudoh[i]=f[i];
#endif
  }

  if(false)
  {
    array2<Complex> f(m,m,align);
    array2<Complex> g(m,m,align);
    array2<Complex> h(m,m,align);
    cconvolution2 convolve(m,f);
    init(f,g);
    seconds();
    convolve.direct(h,f,g);
    sum=seconds();
  
    cout << endl;
    cout << "Direct:" << endl;
    cout << sum-offset/N << endl;
    cout << endl;

    if(m < 5) 
      for(unsigned int i=0; i < m; i++) {
        for(unsigned int j=0; j < m; j++)
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

