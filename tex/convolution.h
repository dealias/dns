#include "fftw++.h"

#ifndef __convolution_h__
#define __convolution_h__ 1

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

class convolution {
protected:
  unsigned int n;
  unsigned int m;
  unsigned int n2;
  unsigned int c;
  rcfft1d *rc;
  crfft1d *cr;
  rcfft1d *rco;
  crfft1d *cro;
  Complex *Zeta;
  Complex Zetam;
  Complex *F;
  Complex *G;
  Complex *B;
  bool even;
public:  
  convolution(unsigned int n, unsigned int m) : n(n), m(m), n2(n/2) {
    rc=new rcfft1d(n);
    cr=new crfft1d(n);
  }
  
  convolution(unsigned int m) : m(m) {
    n=3*m;
    c=m/2;
    
    Zeta=FFTWComplex(c+1);
    
    double arg=2.0*M_PI/n;
    Complex zeta=Complex(cos(arg),sin(arg));
    Complex product=Zeta[0]=1.0;
    for(unsigned int i=1; i <= c; ++i) {
      product *= zeta;
      Zeta[i]=product;
    }
    Zetam=Zeta[c]*Zeta[c];
    even=2.0*c == m;
    if(!even) Zetam *= zeta;
    
    F=FFTWComplex(c+1);
    G=FFTWComplex(c+1);
    B=FFTWComplex(c+1);
    
    rc=new rcfft1d(m,(double *)F,B);
    cr=new crfft1d(m,B,(double *)G);
  }
  
// Need destructor  
  
// Compute H = F (*) G, where F and G are the non-negative Fourier
// components of real functions f and g, respectively. Dealiasing via
// zero-padding is implemented automatically.
//
// Arrays F[n/2+1], g[n/2+1] must be distinct.
// Input F[i] (0 <= i < m), where 3*m <= n, g[i] (0 <= i < n/2).
// Output H[i] = F (*) G  (0 <= i < m), F[i]=f[i], g[i] (0 <= i < n/2).
//
// Array H[n/2+1] can coincide with either F or g, in which case the output H
// subsumes F or g, respectively.

  void fft0(Complex *H, Complex *F, Complex *g) {
    for(unsigned int i=m; i <= n2; i++) F[i]=0.0;
    cr->fft(F);
  
    double ninv=1.0/n;
    for(unsigned int i=0; i < n2; i++)
      H[i]=Complex(F[i].real()*g[i].real()*ninv,F[i].imag()*g[i].imag()*ninv);
	
    rc->fft(H);
  }

// Compute H = F (*) G, where F and G are the non-negative Fourier
// components of real functions f and g, respectively. Dealiasing via
// zero-padding is implemented automatically.
//
// Arrays F[n/2+1], G[n/2+1] must be distinct.
// Input F[i], G[i] (0 <= i < m), where 3*m <= n.
// Output H[i] = F (*) G  (0 <= i < m), F[i]=f[i], G[i]=g[i] (0 <= i < n/2).
//
// Array H[n/2+1] can coincide with either F or G, in which case the output H
// subsumes f or g, respectively.

  void fft(Complex *H, Complex *F, Complex *G) {
    for(unsigned int i=m; i <= n2; i++) G[i]=0.0;
    cr->fft(G);
	
    fft0(H,F,G);
  }	

  // out and in may coincide
  void sym(Complex *in, Complex *out) {
    out[0]=2.0*in[0].real();
    for(int k=1; k <= c; ++k)
      out[k]=in[k]+conj(in[m-k]);
  }

  void unpadded(Complex *h, Complex *f, Complex *g) {
    int stop=even ? c-1 : c;
  
    // r=0:
    double f0=f[0].real();
    double g0=g[0].real();
    
    sym(f,h);

    cr->fft(h,F);
    
    sym(g,h);
    cr->fft(h,G);
    
    for(unsigned int i=0; i < c; i++)
      F[i]=Complex((F[i].real()-f0)*(G[i].real()-g0),
                   (F[i].imag()-f0)*(G[i].imag()-g0));
    F[c]=(F[c].real()-f0)*(G[c].real()-g0);

    rc->fft(F,h);
    
    for(int k=1; k <= stop; ++k)
      h[m-k]=conj(h[k]);

    // r=1:
    B[0]=2.0*f0;
    for(int k=1; k <= c; ++k)
      B[k]=Zeta[k]*(f[k]+conj(Zetam*f[m-k]));

    cr->fft(B,F);
    
    B[0]=2.0*g0;
    for(int k=1; k <= c; ++k)
      B[k]=Zeta[k]*(g[k]+conj(Zetam*g[m-k]));
    cr->fft(B,G);
    
    for(unsigned int i=0; i < c; i++)
      F[i]=Complex((F[i].real()-f0)*(G[i].real()-g0),
                   (F[i].imag()-f0)*(G[i].imag()-g0));
    F[c]=(F[c].real()-f0)*(G[c].real()-g0);

    rc->fft(F,B);
    
    h[0] += B[0].real();
    Complex Zetamc=conj(Zetam);
    for(int k=1; k <= stop; ++k) {
      Complex Bk=multconj(B[k],Zeta[k]);
      h[k] += Bk;
      h[m-k] += multconj(Zetamc,Bk);
    }
    if(even) h[c] += multconj(B[c].real(),Zeta[c]);

    // r=2:
    B[0]=2.0*f0;
    for(int k=1; k <= c; ++k)
      B[k]=multconj(f[k]+Zetam*conj(f[m-k]),Zeta[k]);
    cr->fft(B,F);

    B[0]=2.0*g0;
    for(int k=1; k <= c; ++k)
      B[k]=multconj(g[k]+Zetam*conj(g[m-k]),Zeta[k]);
    cr->fft(B,G);
    
    for(unsigned int i=0; i < c; i++)
      F[i]=Complex((F[i].real()-f0)*(G[i].real()-g0),
                   (F[i].imag()-f0)*(G[i].imag()-g0));
    F[c]=(F[c].real()-f0)*(G[c].real()-g0);

    rc->fft(F,B);
    
    h[0] += B[0].real();
    for(int k=1; k <= stop; ++k) {
      Complex Bk=Zeta[k]*B[k];
      h[k] += Bk;
      h[m-k] += multconj(Zetam,Bk);
    }
    if(even) h[c] += Zeta[c]*B[c].real();

    double ninv=1.0/n;
    for(int i=0; i < m; ++i)
      h[i] *= ninv;
  }
  
// Compute H = F (*) G, where F and G contain the non-negative Fourier
// components of real functions f and g, respectively, via direct convolution
// instead of a Fast Fourier Transform technique.
//
// Input F[i], G[i] (0 <= i < m).
// Output H[i] = F (*) G  (0 <= i < m), F and G unchanged.
//
// Array H[m] must be distinct from F[m] and G[m].

  void direct(Complex *H, Complex *F, Complex *G) {
    for(unsigned int i=0; i < m; i++) {
      Complex sum=0.0;
      for(unsigned int j=0; j <= i; j++) sum += F[j]*G[i-j];
      for(unsigned int j=i+1; j < m; j++) sum += F[j]*conj(G[j-i]);
      for(unsigned int j=1; j < m-i; j++) sum += conj(F[j])*G[i+j];
      H[i]=sum;
    }
  }	

};

#endif
