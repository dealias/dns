// Add namespace to fftw++ and here and cmult-sse2.

#include "fftw++.h"

#include "cmult-sse2.h"

#ifndef __convolution_h__
#define __convolution_h__ 1

extern const double sqrt3;
extern const double hsqrt3;

extern const Complex hSqrt3;
extern const Complex mhsqrt3;
extern const Complex mhalf;
extern const Complex zeta3;

inline unsigned int min(unsigned int a, unsigned int b)
{
  return (a < b) ? a : b;
}

// Build the factored zeta tables.
unsigned int BuildZeta(unsigned int n, unsigned int m,
                       Complex *&ZetaH, Complex *&ZetaL);

// In-place explicitly dealiased complex convolution.
class ExplicitConvolution {
protected:
  unsigned int n,m;
  fft1d *Backwards,*Forwards;
public:  
  
  // u is a temporary array of size n.
  ExplicitConvolution(unsigned int n, unsigned int m, Complex *u) :
    n(n), m(m) {
    Backwards=new fft1d(n,1,u);
    Forwards=new fft1d(n,-1,u);
  }
  
  ~ExplicitConvolution() {
    delete Forwards;
    delete Backwards;
  }    
  
  void padBackwards(Complex *f);
  
  // Compute f (*) g. The distinct input arrays f and g are each of size n 
  // (contents not preserved). The output is returned in f.
  void convolve(Complex *f, Complex *g);
};

// In-place implicitly dealiased complex convolution.
class ImplicitConvolution {
protected:
  unsigned int m;
  unsigned int s;
  fft1d *Backwards,*Forwards;
  fft1d *Backwardso,*Forwardso;
  Complex *ZetaH, *ZetaL;
public:  
  
  // u and v are distinct temporary arrays each of size m.
  ImplicitConvolution(unsigned int m, Complex *u, Complex *v) : m(m) {
    Backwards=new fft1d(m,1,u);
    Forwards=new fft1d(m,-1,u);
    
    Backwardso=new fft1d(m,1,u,v);
    Forwardso=new fft1d(m,-1,u,v);
    
    s=BuildZeta(2*m,m,ZetaH,ZetaL);
  }
  
  ~ImplicitConvolution() {
    deleteAlign(ZetaL);
    deleteAlign(ZetaH);
    delete Forwardso;
    delete Backwardso;
    delete Forwards;
    delete Backwards;
  }
  
  // In-place implicitly dealiased convolution.
  // The input arrays f and g are each of size m (contents not preserved).
  // The output is returned in f.
  // u and v are temporary arrays each of size m.
  void convolve(Complex *f, Complex *g, Complex *u, Complex *v);
};

// Out-of-place direct complex convolution.
class DirectConvolution {
protected:
  unsigned int m;
public:  
  DirectConvolution(unsigned int m) : m(m) {}
  
  void convolve(Complex *h, Complex *f, Complex *g);
};

// In-place explicitly dealiased Hermitian convolution.
class ExplicitHConvolution {
protected:
  unsigned int n,m;
  rcfft1d *rc;
  crfft1d *cr;
public:
  // u is a temporary array of size n.
  ExplicitHConvolution(unsigned int n, unsigned int m, Complex *u) :
    n(n), m(m) {
    rc=new rcfft1d(n,u);
    cr=new crfft1d(n,u);
  }
  
  ~ExplicitHConvolution() {
    delete cr;
    delete rc;
  }
    
  void padBackwards(Complex *f);
  
// Compute f (*) g, where f and g contain the m non-negative Fourier
// components of real functions. Dealiasing is internally implemented via
// explicit zero-padding to size n >= 3*m.
//
// The (distinct) input arrays f and g must each be allocated to size n/2+1
// (contents not preserved). The output is returned in the first m elements
// of f.
  void convolve(Complex *f, Complex *g);
};

// In-place implicitly dealiased Hermitian convolution.
class ImplicitHConvolution {
protected:
  unsigned int m;
  unsigned int c;
  unsigned int s;
  rcfft1d *rc, *rco;
  crfft1d *cr, *cro;
  Complex *ZetaH, *ZetaL;
public:  
  
  // u and v must be each allocated as m/2+1 Complex values.
  ImplicitHConvolution(unsigned int m, Complex *u, Complex *v) : 
    m(m), c(m/2) {

    if(m % 2) {
      std::cerr << "Odd-sized Hermitian convolutions are not implemented;" 
                << std::endl
                << "pad the input data with a final zero element."
                << std::endl;
      _exit(1);
    }
    
    rc=new rcfft1d(m,u);
    cr=new crfft1d(m,u);

    double *U=(double *) u;
    rco=new rcfft1d(m,U,v);
    cro=new crfft1d(m,v,U);
    
    s=BuildZeta(3*m,c,ZetaH,ZetaL);
  }
  
  ~ImplicitHConvolution() {
    deleteAlign(ZetaL);
    deleteAlign(ZetaH);
    delete cro;
    delete rco;
    delete cr;
    delete rc;
  }
  
  // Note: input arrays f and g are destroyed.
  // f and g are each of size m.
  // u and v are temporary arrays each of size m/2+1.
  void convolve(Complex *f, Complex *g, Complex *u, Complex *v);
};
  
// Out-of-place direct Hermitian convolution.
class DirectHConvolution {
protected:
  unsigned int m;
public:  
  DirectHConvolution(unsigned int m) : m(m) {}
  
// Compute h= f (*) g via direct convolution, where f and g contain the m
// non-negative Fourier components of real functions (contents
// preserved). The output of m complex values is returned in the array h,
// which must be distinct from f and g.
  void convolve(Complex *h, Complex *f, Complex *g);
};

// Compute the scrambled virtual m-padded complex Fourier transform of M complex
// vectors, each of length m.
// The arrays in and out (which may coincide), along with the array u, must
// be allocated as Complex[M*m].
//
//   fftpad fft(m,M,stride);
//   fft.backwards(in,u);
//   fft.forwards(in,u);
//
// Notes:
//   stride is the spacing between the elements of each Complex vector.
//
class fftpad {
  unsigned int m;
  unsigned int M;
  unsigned int stride;
  unsigned int dist;
  unsigned int s;
  mfft1d *Backwards;
  mfft1d *Forwards;
  Complex *ZetaH, *ZetaL;
public:  
  fftpad(unsigned int m, unsigned int M,
         unsigned int stride, Complex *f) : m(m), M(M), stride(stride) {
    Backwards=new mfft1d(m,1,M,stride,1,f);
    Forwards=new mfft1d(m,-1,M,stride,1,f);
    
    s=BuildZeta(2*m,m,ZetaH,ZetaL);
  }
  
  ~fftpad() {
    deleteAlign(ZetaL);
    deleteAlign(ZetaH);
    delete Forwards;
    delete Backwards;
  }
  
  void backwards(Complex *f, Complex *u);
  void forwards(Complex *f, Complex *u);
};
  
// Compute the scrambled virtual m-padded complex Fourier transform of M complex
// vectors, each of length 2m-1 with the origin at index m-1
// (i.e. physical wavenumber k=-m+1 to k=m-1).
// The arrays in and out (which may coincide) must be allocated as
// Complex[M*(2m-1)]. The array u must be allocated as Complex[M*(m+1)].
//
//   fft0pad fft(m,M,stride);
//   fft.backwards(in,u);
//   fft.forwards(in,u);
//
// Notes:
//   stride is the spacing between the elements of each Complex vector.
//
class fft0pad {
  unsigned int m;
  unsigned int M;
  unsigned int s;
  unsigned int stride;
  mfft1d *Backwards;
  mfft1d *Forwards;
  Complex *ZetaH, *ZetaL;
public:  
  fft0pad(unsigned int m, unsigned int M, unsigned int stride, Complex *u)
    : m(m), M(M), stride(stride) {
    Backwards=new mfft1d(m,1,M,stride,1,u);
    Forwards=new mfft1d(m,-1,M,stride,1,u);
    
    s=BuildZeta(3*m,m,ZetaH,ZetaL);
  }
  
  ~fft0pad() {
    deleteAlign(ZetaL);
    deleteAlign(ZetaH);
    delete Forwards;
    delete Backwards;
  }
  
  void backwards(Complex *f, Complex *u);
  void forwards(Complex *f, Complex *u);
};
  
// In-place explicitly dealiased 2D complex convolution.
class ExplicitConvolution2 {
protected:
  unsigned int nx,ny;
  unsigned int mx,my;
  bool prune; // Skip Fourier transforming rows containing all zeroes?
  mfft1d *xBackwards, *xForwards;
  mfft1d *yBackwards, *yForwards;
  fft2d *Backwards, *Forwards;
public:
  ExplicitConvolution2(unsigned int nx, unsigned int ny,
                       unsigned int mx, unsigned int my,
                       Complex *f, bool prune=false) :
    nx(nx), ny(ny), mx(mx), my(my), prune(prune) {
    if(prune) {
      xBackwards=new mfft1d(nx,1,my,ny,1,f);
      yBackwards=new mfft1d(ny,1,nx,1,ny,f);
      yForwards=new mfft1d(ny,-1,nx,1,ny,f);
      xForwards=new mfft1d(nx,-1,my,ny,1,f);
    } else {
      Backwards=new fft2d(nx,ny,1,f);
      Forwards=new fft2d(nx,ny,-1,f);
    }
  }
  
  ~ExplicitConvolution2() {
    if(prune) {
      delete xForwards;
      delete yForwards;
      delete yBackwards;
      delete xBackwards;
    } else {
      delete Forwards;
      delete Backwards;
    }
  }    
  
  void padBackwards(Complex *f);
  void convolve(Complex *f, Complex *g);
};

// In-place implicitly dealiased 2D complex convolution.
class ImplicitConvolution2 {
protected:
  unsigned int mx,my;
  fftpad *xfftpad;
  ImplicitConvolution *yconvolve;
public:  
  // u1 and v1 are temporary arrays of size my.
  // u2 is a temporary array of size mx*my.
  ImplicitConvolution2(unsigned int mx, unsigned int my,
                       Complex *u1, Complex *v1, Complex *u2) : mx(mx), my(my) {
    xfftpad=new fftpad(mx,my,my,u2);
    yconvolve=new ImplicitConvolution(my,u1,v1);
  }
  
  ~ImplicitConvolution2() {
    delete yconvolve;
    delete xfftpad;
  }
  
  // The distinct input arrays f and g must be allocated as Complex[mx*my]. 
  // u1 and v1 are temporary arrays of size my.
  // u2 and v2 are temporary arrays of size mx*my.
  // The output is returned in f.
  void convolve(Complex *f, Complex *g, Complex *u1, Complex *v1,
                Complex *u2, Complex *v2) {
    xfftpad->backwards(f,u2);
    xfftpad->backwards(g,v2);

    unsigned int mxy=mx*my;
    for(unsigned int i=0; i < mxy; i += my)
      yconvolve->convolve(f+i,g+i,u1,v1);
    for(unsigned int i=0; i < mxy; i += my)
      yconvolve->convolve(u2+i,v2+i,u1,v1);
    
    xfftpad->forwards(f,u2);
  }
};

// Out-of-place direct 2D complex convolution.
class DirectConvolution2 {
protected:  
  unsigned int mx,my;
public:
  DirectConvolution2(unsigned int mx, unsigned int my) : mx(mx), my(my) {}
  
  void convolve(Complex *h, Complex *f, Complex *g);
};

// Enforce Hermiticity by symmetrizing 2D data on the X axis.
inline void HermitianSymmetrizeX(unsigned int mx, unsigned int my,
                                 unsigned int xorigin, Complex *f)
{
  unsigned int offset=xorigin*my;
  unsigned int stop=mx*my;
  for(unsigned int i=my; i < stop; i += my)
    f[offset-i]=conj(f[offset+i]);
}

// Enforce Hermiticity by symmetrizing D data on the X axis.
inline void HermitianSymmetrizeX(unsigned int mx, unsigned int my,
                                 unsigned int mz, unsigned int ny,
                                 unsigned int xorigin, unsigned int yorigin, 
                                 Complex *f)
{
  for(unsigned int i=1; i < mx; ++i)
    f[((xorigin-i)*ny+yorigin)*mz]=conj(f[((xorigin+i)*ny+yorigin)*mz]);
}

// Enforce Hermiticity by symmetrizing D data on the XY plane.
inline void HermitianSymmetrizeXY(unsigned int mx, unsigned int my,
                                  unsigned int mz, unsigned int ny,
                                  unsigned int xorigin, unsigned int yorigin, 
                                  Complex *f)
{
  for(unsigned int i=1; i < mx; ++i)
    f[((xorigin-i)*ny+yorigin)*mz]=conj(f[((xorigin+i)*ny+yorigin)*mz]);
  
  for(unsigned int i=0; i < mx; ++i)
    for(unsigned int j=1; j < my; ++j) {
      f[((xorigin-i)*ny+yorigin-j)*mz]=conj(f[((xorigin+i)*ny+(yorigin+j))*mz]);
      f[((xorigin+i)*ny+yorigin-j)*mz]=conj(f[((xorigin-i)*ny+(yorigin+j))*mz]);
    }
}

// In-place explicitly dealiased 2D Hermitian convolution.
class ExplicitHConvolution2 {
protected:
  unsigned int nx,ny;
  unsigned int mx,my;
  bool prune; // Skip Fourier transforming rows containing all zeroes?
  mfft1d *xBackwards;
  mfft1d *xForwards;
  mcrfft1d *yBackwards;
  mrcfft1d *yForwards;
  crfft2d *Backwards;
  rcfft2d *Forwards;
public:  
  ExplicitHConvolution2(unsigned int nx, unsigned int ny, 
                        unsigned int mx, unsigned int my, Complex *f,
                        bool pruned=false) :
    nx(nx), ny(ny), mx(mx), my(my), prune(pruned) {
    unsigned int nyp=ny/2+1;
    // Odd nx requires interleaving of shift with x and y transforms.
    unsigned int My=my;
    if(nx % 2) {
      if(!prune) My=nyp;
      prune=true;
    }

    if(prune) {
      xBackwards=new mfft1d(nx,1,My,nyp,1,f);
      yBackwards=new mcrfft1d(ny,nx,1,nyp,f);
      yForwards=new mrcfft1d(ny,nx,1,nyp,f);
      xForwards=new mfft1d(nx,-1,My,nyp,1,f);
    } else {
      Backwards=new crfft2d(nx,ny,f);
      Forwards=new rcfft2d(nx,ny,f);
    }
  }
  
  ~ExplicitHConvolution2() {
    if(prune) {
      delete xForwards;
      delete yForwards;
      delete yBackwards;
      delete xBackwards;
    } else {
      delete Forwards;
      delete Backwards;
    }
  }    
  
  void padBackwards(Complex *f);
  void convolve(Complex *f, Complex *g);
};

// In-place implicitly dealiased 2D Hermitian convolution.
class ImplicitHConvolution2 {
protected:
  unsigned int mx,my;
  fft0pad *xfftpad;
  ImplicitHConvolution *yconvolve;
public:  
  // u1 and v1 are temporary arrays of size my/2+1.
  // u2 is a temporary array of size (mx+1)*my;
  ImplicitHConvolution2(unsigned int mx, unsigned int my,
                        Complex *u1, Complex *v1, Complex *u2) :
    mx(mx), my(my) {
    xfftpad=new fft0pad(mx,my,my,u2);
    yconvolve=new ImplicitHConvolution(my,u1,v1);
  }
  
  ~ImplicitHConvolution2() {
    delete yconvolve;
    delete xfftpad;
  }
  
  // The distinct input arrays f and g must be allocated as Complex[(2mx-1)*my] 
  // (not preserved).
  // u1 and v1 are temporary arrays allocated as Complex[my/2+1].
  // u2 and v2 are temporary arrays allocated as Complex [(mx+1)*my];
  // The output is returned in f.
  void convolve(Complex *f, Complex *g, Complex *u1, Complex *v1,
                Complex *u2, Complex *v2) {
    unsigned int xorigin=mx-1;
    
    HermitianSymmetrizeX(mx,my,xorigin,f);
    HermitianSymmetrizeX(mx,my,xorigin,g);
    
    xfftpad->backwards(f,u2);
    xfftpad->backwards(g,v2);

    unsigned int mf=(2*mx-1)*my;
    for(unsigned int i=0; i < mf; i += my)
      yconvolve->convolve(f+i,g+i,u1,v1);
    unsigned int mu=(mx+1)*my;
    for(unsigned int i=0; i < mu; i += my)
      yconvolve->convolve(u2+i,v2+i,u1,v1);
    
    xfftpad->forwards(f,u2);
  }
};

// Out-of-place direct 2D Hermitian convolution.
class DirectHConvolution2 {
protected:  
  unsigned int mx,my;
public:
  DirectHConvolution2(unsigned int mx, unsigned int my) : mx(mx), my(my) {}
  
  void convolve(Complex *h, Complex *f, Complex *g);
};

// In-place explicitly dealiased 3D complex convolution.
class ExplicitConvolution3 {
protected:
  unsigned int nx,ny,nz;
  unsigned int mx,my,mz;
  bool prune; // Skip Fourier transforming rows containing all zeroes?
  mfft1d *xBackwards, *xForwards;
  mfft1d *yBackwards, *yForwards;
  mfft1d *zBackwards, *zForwards;
  fft3d *Backwards, *Forwards;
public:
  ExplicitConvolution3(unsigned int nx, unsigned int ny, unsigned int nz,
                       unsigned int mx, unsigned int my, unsigned int mz,
                       Complex *f, bool prune=false) :
    nx(nx), ny(ny), nz(nz), mx(mx), my(my), mz(mz), prune(prune) {
    unsigned int nxy=nx*ny;
    unsigned int nyz=ny*nz;
    if(prune) {
      xBackwards=new mfft1d(nx,1,nyz,nyz,1,f);
      yBackwards=new mfft1d(ny,1,mz,nz,1,f);
      zBackwards=new mfft1d(nz,1,nxy,1,nz,f);
      zForwards=new mfft1d(nz,-1,nxy,1,nz,f);
      yForwards=new mfft1d(ny,-1,mz,nz,1,f);
      xForwards=new mfft1d(nx,-1,nyz,nyz,1,f);
    } else {
      Backwards=new fft3d(nx,ny,nz,1,f);
      Forwards=new fft3d(nx,ny,nz,-1,f);
    }
  }
  
  ~ExplicitConvolution3() {
    if(prune) {
      delete xForwards;
      delete yForwards;
      delete zForwards;
      delete zBackwards;
      delete yBackwards;
      delete xBackwards;
    } else {
      delete Forwards;
      delete Backwards;
    }
  }    
  
  void padBackwards(Complex *f);
  void convolve(Complex *f, Complex *g);
};

// In-place implicitly dealiased 3D complex convolution.
class ImplicitConvolution3 {
protected:
  unsigned int mx,my,mz;
  fftpad *xfftpad;
  ImplicitConvolution2 *yzconvolve;
public:  
  // u1 and v1 are temporary arrays of size mz.
  // u2 is a temporary array of size my*mz.
  // u3 is a temporary array of size mx*my*mz.
  ImplicitConvolution3(unsigned int mx, unsigned int my, unsigned int mz,
                       Complex *u1, Complex *v1, Complex *u2, Complex *u3) : 
    mx(mx), my(my), mz(mz) {
    xfftpad=new fftpad(mx,my*mz,my*mz,u3);
    yzconvolve=new ImplicitConvolution2(my,mz,u1,v1,u2);
  }
  
  ~ImplicitConvolution3() {
    delete yzconvolve;
    delete xfftpad;
  }
  
  // The distinct input arrays f and g must be allocated as Complex[mx*my*mz]. 
  // u1 and v1 are temporary arrays of size mz.
  // u2 and v2 are temporary arrays of size my*mz.
  // u3 and v3 are temporary arrays of size mx*my*mz.
  // The output is returned in f.
  void convolve(Complex *f, Complex *g, Complex *u1, Complex *v1,
                Complex *u2, Complex *v2, Complex *u3, Complex *v3) {
    xfftpad->backwards(f,u3);
    xfftpad->backwards(g,v3);

    unsigned int myz=my*mz;
    unsigned int mxyz=mx*myz;
    for(unsigned int i=0; i < mxyz; i += myz)
      yzconvolve->convolve(f+i,g+i,u1,v1,u2,v2);
    for(unsigned int i=0; i < mxyz; i += myz)
      yzconvolve->convolve(u3+i,v3+i,u1,v1,u2,v2);
    
    xfftpad->forwards(f,u3);
  }
};

// Out-of-place direct 3D complex convolution.
class DirectConvolution3 {
protected:  
  unsigned int mx,my,mz;
  unsigned int myz;
public:
  DirectConvolution3(unsigned int mx, unsigned int my, unsigned int mz) : 
    mx(mx), my(my), mz(mz), myz(my*mz) {}
  
  void convolve(Complex *h, Complex *f, Complex *g);
};

// In-place implicitly dealiased 3D Hermitian convolution.
class ImplicitHConvolution3 {
protected:
  unsigned int mx,my,mz;
  fft0pad *xfftpad;
  ImplicitHConvolution2 *yzconvolve;
public:  
  // u1 and v1 are temporary arrays of size mz/2+1.
  // u2 is a temporary array of size (my+1)*mz.
  // u3 is a temporary array of size (mx+1)*(2my-1)*mz.
  ImplicitHConvolution3(unsigned int mx, unsigned int my, unsigned int mz,
                        Complex *u1, Complex *v1, Complex *u2, Complex *u3) : 
    mx(mx), my(my), mz(mz) {
    unsigned int nymz=(2*my-1)*mz;
    xfftpad=new fft0pad(mx,nymz,nymz,u3);
    yzconvolve=new ImplicitHConvolution2(my,mz,u1,v1,u2);
  }
  
  ~ImplicitHConvolution3() {
    delete yzconvolve;
    delete xfftpad;
  }
  
  // The distinct input arrays f and g must be allocated as 
  // Complex[(2mx-1)*(2my-1)*mz] (not preserved). 
  // u1 and v1 are temporary arrays allocated as Complex[mz/2+1].
  // u2 and v2 are temporary arrays allocated as Complex[(my+1)*mz].
  // u3 and v3 are temporary arrays allocated as Complex[(mx+1)*(2my-1)*mz].
  // The output is returned in f.
  void convolve(Complex *f, Complex *g, Complex *u1, Complex *v1,
                Complex *u2, Complex *v2, Complex *u3, Complex *v3) {
    unsigned int xorigin=mx-1;
    unsigned int yorigin=my-1;
    unsigned int ny=2*my-1;
    
    HermitianSymmetrizeX(mx,my,mz,ny,xorigin,yorigin,f);
    HermitianSymmetrizeX(mx,my,mz,ny,xorigin,yorigin,g);
    
    xfftpad->backwards(f,u3);
    xfftpad->backwards(g,v3);

    unsigned int incr=ny*mz;
    unsigned int fstop=(2*mx-1)*incr;
    for(unsigned int i=0; i < fstop; i += incr)
      yzconvolve->convolve(f+i,g+i,u1,v1,u2,v2);
    unsigned int ustop=(mx+1)*incr;
    for(unsigned int i=0; i < ustop; i += incr)
      yzconvolve->convolve(u3+i,v3+i,u1,v1,u2,v2);
    
    xfftpad->forwards(f,u3);
  }
};

// Out-of-place direct 3D Hermitian convolution.
class DirectHConvolution3 {
protected:  
  unsigned int mx,my,mz;
public:
  DirectHConvolution3(unsigned int mx, unsigned int my, unsigned int mz) : 
    mx(mx), my(my), mz(mz) {}
  
  void convolve(Complex *h, Complex *f, Complex *g);
};

// In-place explicitly dealiased Hermitian biconvolution.
class ExplicitHBiConvolution {
protected:
  unsigned int n;
  unsigned int m;
  rcfft1d *rc;
  crfft1d *cr;
public:
  // u is a temporary array of size n.
  ExplicitHBiConvolution(unsigned int n, unsigned int m, Complex *u) :
    n(n), m(m) {
    rc=new rcfft1d(n,u);
    cr=new crfft1d(n,u);
  }
  
  ~ExplicitHBiConvolution() {
    delete cr;
    delete rc;
  }
    
  void padBackwards(Complex *f);
  
// Compute the biconvolution of f, and g, and h, where f, and g, and h
// contain the m non-negative Fourier components of real
// functions. Dealiasing is internally implemented via explicit
// zero-padding to size n >= 3*m. The (distinct) input arrays f, g, and h
// must each be allocated to size n/2+1 (contents not preserved).
// The output is returned in the first m elements of f.
  void convolve(Complex *f, Complex *g, Complex *h);
};

// In-place implicitly dealiased Hermitian biconvolution.
class ImplicitHBiConvolution {
protected:
  unsigned int m;
  unsigned int s;
  rcfft1d *rc, *rco;
  crfft1d *cr, *cro;
  Complex *ZetaH, *ZetaL;
public:  
  
  // u and v are distinct temporary arrays each of size m+1.
  ImplicitHBiConvolution(unsigned int m, Complex *u, Complex *v) : 
    m(m) {
    unsigned int twom=2*m;
    
    rc=new rcfft1d(twom,u);
    cr=new crfft1d(twom,u);
    
    double *U=(double *) u;
    rco=new rcfft1d(twom,U,v);
    cro=new crfft1d(twom,v,U);
    
    s=BuildZeta(4*m,m,ZetaH,ZetaL);
  }
  
  ~ImplicitHBiConvolution() {
    deleteAlign(ZetaL);
    deleteAlign(ZetaH);
    delete cro;
    delete rco;
    delete cr;
    delete rc;
  }
  
  // In-place implicitly dealiased convolution.
  // The input arrays f, g, and h are each of size m+1 (contents not preserved).
  // The output is returned in f.
  // u, v, and w are temporary arrays each of size m+1.
  void convolve(Complex *f, Complex *g, Complex *h,
                Complex *u, Complex *v, Complex *w);
};

// Out-of-place direct 1D Hermitian biconvolution.
class DirectHBiConvolution {
protected:  
  unsigned int m;
public:
  DirectHBiConvolution(unsigned int m) : m(m) {}
  
  void convolve(Complex *h, Complex *e, Complex *f, Complex *g);
};

// Compute the scrambled virtual 2m-padded complex Fourier transform of M
// complex vectors, each of length 2m with the Fourier origin at index m.
// The arrays in and out (which may coincide), along
// with the array u, must be allocated as Complex[M*2m].
//
//   fft0bipad fft(m,M,stride);
//   fft.backwards(in,u);
//   fft.forwards(in,u);
//
// Notes:
//   stride is the spacing between the elements of each Complex vector.
//
class fft0bipad {
  unsigned int m;
  unsigned int M;
  unsigned int stride;
  unsigned int s;
  mfft1d *Backwards;
  mfft1d *Forwards;
  Complex *ZetaH, *ZetaL;
public:  
  fft0bipad(unsigned int m, unsigned int M, unsigned int stride,
            Complex *f) : m(m), M(M), stride(stride) {
    unsigned int twom=2*m;
    Backwards=new mfft1d(twom,1,M,stride,1,f);
    Forwards=new mfft1d(twom,-1,M,stride,1,f);
    
    s=BuildZeta(4*m,twom,ZetaH,ZetaL);
  }
  
  ~fft0bipad() {
    deleteAlign(ZetaL);
    deleteAlign(ZetaH);
    delete Forwards;
    delete Backwards;
  }
  
  void backwards(Complex *f, Complex *u);
  void forwards(Complex *f, Complex *u);
};
  
// In-place explicitly dealiased 2D Hermitian biconvolution.
class ExplicitHBiConvolution2 {
protected:
  unsigned int nx,ny;
  unsigned int mx,my;
  bool prune; // Skip Fourier transforming rows containing all zeroes?
  mfft1d *xBackwards;
  mfft1d *xForwards;
  mcrfft1d *yBackwards;
  mrcfft1d *yForwards;
  crfft2d *Backwards;
  rcfft2d *Forwards;
public:  
  ExplicitHBiConvolution2(unsigned int nx, unsigned int ny, 
                          unsigned int mx, unsigned int my, Complex *f,
                          bool pruned=false) :
    nx(nx), ny(ny), mx(mx), my(my), prune(pruned) {
    unsigned int nyp=ny/2+1;
    // Odd nx requires interleaving of shift with x and y transforms.
    unsigned int My=my;
    if(nx % 2) {
      if(!prune) My=nyp;
      prune=true;
    }
    
    if(prune) {
      xBackwards=new mfft1d(nx,1,My,nyp,1,f);
      yBackwards=new mcrfft1d(ny,nx,1,nyp,f);
      yForwards=new mrcfft1d(ny,nx,1,nyp,f);
      xForwards=new mfft1d(nx,-1,My,nyp,1,f);
    } else {
      Backwards=new crfft2d(nx,ny,f);
      Forwards=new rcfft2d(nx,ny,f);
    }
  }
  
  ~ExplicitHBiConvolution2() {
    if(prune) {
      delete xForwards;
      delete yForwards;
      delete yBackwards;
      delete xBackwards;
    } else {
      delete Forwards;
      delete Backwards;
    }
  }    
  
  void padBackwards(Complex *f);
  void convolve(Complex *f, Complex *g, Complex *h);
};

// In-place implicitly dealiased 2D Hermitian biconvolution.
class ImplicitHBiConvolution2 {
protected:
  unsigned int mx,my;
  fft0bipad *xfftpad;
  ImplicitHBiConvolution *yconvolve;
public:  
  // u1 and v1 are temporary arrays of size my+1.
  // u2 is a temporary array of size 2mx*(my+1);
  ImplicitHBiConvolution2(unsigned int mx, unsigned int my,
                          Complex *u1, Complex *v1, Complex *u2) :
    mx(mx), my(my) {
    xfftpad=new fft0bipad(mx,my,my+1,u2);
    yconvolve=new ImplicitHBiConvolution(my,u1,v1);
  }
  
  ~ImplicitHBiConvolution2() {
    delete yconvolve;
    delete xfftpad;
  }
  
  // The distinct input arrays f, g, and h must be allocated as 
  // Complex[2mx*(my+1)] (not preserved).
  // u1, v1, and w1 are temporary arrays allocated as Complex[my+1].
  // u2, v2, and w2 are temporary arrays allocated as Complex [2mx*(my+1)];
  // The output is returned in f.
  void convolve(Complex *f, Complex *g, Complex *h,
                Complex *u1, Complex *v1, Complex *w1,
                Complex *u2, Complex *v2, Complex *w2) {
    unsigned int my1=my+1;
    HermitianSymmetrizeX(mx,my1,mx,f);
    HermitianSymmetrizeX(mx,my1,mx,g);
    HermitianSymmetrizeX(mx,my1,mx,h);
    
    xfftpad->backwards(f,u2);
    xfftpad->backwards(g,v2);
    xfftpad->backwards(h,w2);

    unsigned int stop=2*mx*my1;

    for(unsigned int i=0; i < stop; i += my1)
      yconvolve->convolve(f+i,g+i,h+i,u1,v1,w1);
    for(unsigned int i=0; i < stop; i += my1)
      yconvolve->convolve(u2+i,v2+i,w2+i,u1,v1,w1);
    
    xfftpad->forwards(f,u2);
  }
};

// Out-of-place direct 2D Hermitian biconvolution.
class DirectHBiConvolution2 {
protected:  
  unsigned int mx,my;
public:
  DirectHBiConvolution2(unsigned int mx, unsigned int my) : mx(mx), my(my) {}
  
  void convolve(Complex *h, Complex *e, Complex *f, Complex *g);
};

#endif
