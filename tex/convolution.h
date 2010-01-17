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
  Complex *F;
  Complex *G;
  Complex *B;
public:  
  convolution(unsigned int n, unsigned int m) : n(n), m(m), n2(n/2) {
    rc=new rcfft1d(n);
    cr=new crfft1d(n);
  }
  
  convolution(unsigned int m) : m(m) {
    n=3*m;
    c=m/2;
    
    Zeta=FFTWComplex(m+1);
    Complex product=1.0;
    if(n > 0) {
      Zeta[0]=product;
      double arg=2.0*M_PI/n;
      Complex zeta=Complex(cos(arg),sin(arg));
      for(unsigned int i=1; i <= m; ++i) {
        product *= zeta;
        Zeta[i]=product;
      }
    }
    
    F=FFTWComplex(m);
    G=FFTWComplex(m);
    B=FFTWComplex(m);
    
    rc=new rcfft1d(m,(double *)F,G);
    cr=new crfft1d(m,F,(double *)G);
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
    bool even=2.0*c == m;
    int stop=c;
    if(even) --stop;
  
    // r=0:
    double F0=f[0].real();
    double G0=g[0].real();
    
    sym(f,h);

    cr->fft(h,F);
    
    sym(g,h);
    cr->fft(h,G);
    
    for(unsigned int i=0; i < c; i++)
      F[i]=Complex((F[i].real()-F0)*(G[i].real()-G0),
                   (F[i].imag()-F0)*(G[i].imag()-G0));
    F[c]=(F[c].real()-F0)*(G[c].real()-G0);

    rc->fft(F,h);
    
    for(int k=1; k <= stop; ++k)
      h[m-k]=conj(h[k]);

    // r=1:
    F[0]=f[0];
    G[0]=g[0];
    for(int k=1; k < m; ++k) {
      Complex Zetak=Zeta[k];
      F[k]=Zetak*f[k];
      G[k]=Zetak*g[k];
    }
  
    F0=F[0].real();
    G0=G[0].real();
    
    sym(F,B);
    cr->fft(B,F);
    
    sym(G,B);
    cr->fft(B,G);
    
    for(unsigned int i=0; i < c; i++)
      F[i]=Complex((F[i].real()-F0)*(G[i].real()-G0),
                   (F[i].imag()-F0)*(G[i].imag()-G0));
    F[c]=(F[c].real()-F0)*(G[c].real()-G0);

    rc->fft(F,B);
    
    h[0] += B[0].real();
    Complex Zetamc=conj(Zeta[m]);
    for(int k=1; k <= stop; ++k) {
      Complex Bk=multconj(B[k],Zeta[k]);
      h[k] += Bk;
      h[m-k] += multconj(Zetamc,Bk);
    }
    if(even) h[c] += multconj(B[c].real(),Zeta[c]);

    // r=2:
    F[0]=f[0];
    G[0]=g[0];
    for(int k=1; k < m; ++k) {
      Complex Zetamk=Zeta[k];
      F[k]=multconj(f[k],Zetamk);
      G[k]=multconj(g[k],Zetamk);
    }

    F0=F[0].real();
    G0=G[0].real();
    
    sym(F,B);
    cr->fft(B,F);
    
    sym(G,B);
    cr->fft(B,G);
    
    for(unsigned int i=0; i < c; i++)
      F[i]=Complex((F[i].real()-F0)*(G[i].real()-G0),
                   (F[i].imag()-F0)*(G[i].imag()-G0));
    F[c]=(F[c].real()-F0)*(G[c].real()-G0);

    rc->fft(F,B);
    
    h[0] += B[0].real();
    Complex Zetam=Zeta[m];
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
