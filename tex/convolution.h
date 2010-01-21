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
  unsigned int c;
  rcfft1d *rc;
  crfft1d *cr;
  double *F,*G;
  Complex *A,*B;
  Complex zeta;
public:  
  convolution(unsigned int n, unsigned int m, bool compact=false) :
    n(n), m(m) {
    if(compact) {
      F=NULL;
      G=NULL;
      rc=new rcfft1d(n);
      cr=new crfft1d(n);
    } else {
      F=FFTWdouble(n);
      G=FFTWdouble(n+2);
      rc=new rcfft1d(n,F,(Complex *) G);
      cr=new crfft1d(n,(Complex *) G,F);
    }
  }
  
  convolution(unsigned int m) : m(m) {
    n=3*m;
    c=m/2;
    
    double arg=2.0*M_PI/n;
    zeta=Complex(cos(arg),sin(arg));

    // Work arrays:
    A=FFTWComplex(c+1);
    B=FFTWComplex(c+1);
    double *b=(double *) B;
    
    rc=new rcfft1d(m,b,A);
    cr=new crfft1d(m,A,b);
  }
  
// Need destructor  
  
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

  void fft(Complex *h, Complex *f, Complex *g) {
    unsigned int n2=n/2;
    double ninv=1.0/n;
    if(F) {
      for(unsigned int i=m; i <= n2; i++) f[i]=0.0;
      cr->fft(f,F);
  
      for(unsigned int i=m; i <= n2; i++) g[i]=0.0;
      cr->fft(g,G);
      
      for(unsigned int i=0; i < n; ++i)
        F[i] *= G[i]*ninv;
	
      rc->fft(F,h);
    } else {
      for(unsigned int i=m; i <= n2; i++) f[i]=0.0;
      cr->fft(f);
  
      for(unsigned int i=m; i <= n2; i++) g[i]=0.0;
      cr->fft(g);
	
      double *F=(double *) f;
      double *G=(double *) g;
      double *H=(double *) h;
      for(unsigned int i=0; i < n; ++i)
        H[i]=F[i]*G[i]*ninv;
	
      rc->fft(h);
    }
  }
  
  // Note: input arrays f and g are destroyed.
  void unpadded(Complex *h, Complex *f, Complex *g) {
    double f0=f[0].real();
    double g0=g[0].real();

    h[0]=A[0]=f0;
    B[0]=g0;
    Complex *C=h+c-1;
    Complex C1=g[1]+conj(g[m-1]);
    static const Complex Zetamc=Complex(-0.5,-0.5*sqrt(3.0));
    Complex Zetak=zeta;
    for(int k=1; k <= c; ++k) {
      Complex fk=f[k];
      Complex fmk=conj(f[m-k]);
      h[k]=fk+fmk;
      const Complex Zetakm=Zetak*Zetamc;
      f[k]=fk*Zetak+fmk*Zetakm;
      A[k]=multconj(fk,Zetak)+multconj(fmk,Zetakm);
      Complex gk=g[k];
      Complex gmk=conj(g[m-k]);
      C[k]=gk+gmk;
      g[k]=gk*Zetak+gmk*Zetakm;
      B[k]=multconj(gk,Zetak)+multconj(gmk,Zetakm);
      Zetak *= zeta;
    }
    
    Complex fc=f[c];
    double *F=(double *)(f+c);
    double *G=(double *) A;
    
    // r=-1:
    cr->fft(A,F);
    cr->fft(B,G);
    for(unsigned int i=0; i < m; i++)
      F[i] *= G[i];
    rc->fft(F,B);
    
    // r=0:
    cr->fft(h,F);
    C[0]=g0;
    C[1]=C1;
    cr->fft(C,G);
    for(unsigned int i=0; i < m; i++)
      F[i] *= G[i];
    rc->fft(F,h);
    
    // r=1:
    f[c]=fc;
    cr->fft(f,F);
    cr->fft(g,G);
    for(unsigned int i=0; i < m; i++)
      F[i] *= G[i];
    rc->fft(F,g);
    
    double ninv=1.0/n;
    h[0]=(h[0].real()+g[0].real()+B[0].real())*ninv;
    Zetak=zeta*ninv;
    int stop=m-c-1;
    for(int k=1; k <= stop; ++k) {
      Complex gk=multconj(g[k],Zetak);
      Complex Bk=Zetak*B[k];
      Zetak *= zeta;
      h[m-k]=conj(h[k]*ninv+multconj(gk,Zetamc)+Zetamc*Bk);
      h[k]=h[k]*ninv+gk+Bk;
    }
    if(2*c == m) 
      h[c]=h[c].real()*ninv+g[c].real()*conj(Zetak)+B[c].real()*Zetak;
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
