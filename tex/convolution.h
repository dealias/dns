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
  double *F,*G;
  Complex *B,*C;
  Complex zeta;
  bool even;
public:  
  convolution(unsigned int n, unsigned int m) : n(n), m(m), n2(n/2) {
    F=FFTWdouble(n);
    G=FFTWdouble(n+2);
    rc=new rcfft1d(n,F,(Complex *) G);
    cr=new crfft1d(n,(Complex *) G,F);
  }
  
  convolution(unsigned int m) : m(m) {
    n=3*m;
    c=m/2;
    
    double arg=2.0*M_PI/n;
    zeta=Complex(cos(arg),sin(arg));

    even=2.0*c == m;
    
    // Work arrays:
    F=FFTWdouble(m);
    B=FFTWComplex(c+1);
    C=FFTWComplex(c+1);
    
    rc=new rcfft1d(m,F,B);
    cr=new crfft1d(m,B,F);
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
    for(unsigned int i=m; i <= n2; i++) f[i]=0.0;
    cr->fft(f,F);
  
    for(unsigned int i=m; i <= n2; i++) g[i]=0.0;
    cr->fft(g,G);
	
    double ninv=1.0/n;
    for(unsigned int i=0; i < n; ++i)
      F[i] *= G[i]*ninv;
	
    rc->fft(F,h);
  }
  
  void unpadded(Complex *h, Complex *f, Complex *g) {
    int stop=even ? c-1 : c;
    
    Complex f0(f[0].real(),-f[0].real());
    Complex g0(g[0].real(),-g[0].real());

    Complex *A=h+c-1;
    
    h[0]=f[0]=g[0]=f0;
    B[0]=C[0]=g0;
    static const Complex Zetamc=Complex(-0.5,-0.5*sqrt(3.0));
    Complex A1=g[1]+conj(g[m-1]);
    
    Complex Zetak=zeta;
    for(int k=1; k <= c; ++k) {
      Complex fk=f[k];
      Complex fmk=conj(f[m-k]);
      h[k]=fk+fmk;
      const Complex Zetakm=Zetak*Zetamc;
      f[k]=fk*Zetak+fmk*Zetakm;
      Complex gk=g[k];
      Complex gmk=conj(g[m-k]);
      g[k]=multconj(fk,Zetak)+multconj(fmk,Zetakm);
      A[k]=gk+gmk;
      B[k]=gk*Zetak+gmk*Zetakm;
      C[k]=multconj(gk,Zetak)+multconj(gmk,Zetakm);
      Zetak *= zeta;
    }
    
    Complex fcm1=f[c-1];
    Complex fc=f[c];
    double *H=(double *)(f+c-1);
    
    // r=0:
    cr->fft(h,F);
    A[0]=g0;
    A[1]=A1;
    cr->fft(A,H);
    for(unsigned int i=0; i < m; i++)
      F[i] *= H[i];
    rc->fft(F,h);
    
    // r=1:
    f[c-1]=fcm1;
    f[c]=fc;
    cr->fft(f,F);
    cr->fft(B,H);
    for(unsigned int i=0; i < m; i++)
      F[i] *= H[i];
    rc->fft(F,B);
    
    // r=2:
    cr->fft(g,F);
    cr->fft(C,H);
    for(unsigned int i=0; i < m; i++)
      F[i] *= H[i];
    rc->fft(F,C);
    
    double ninv=1.0/n;
    h[0]=(h[0].real()+B[0].real()+C[0].real())*ninv;
    Zetak=zeta*ninv;
    for(int k=1; k <= stop; ++k) {
      Complex Bk=multconj(B[k],Zetak);
      Complex Ck=Zetak*C[k];
      Zetak *= zeta;
      h[m-k]=conj(h[k]*ninv+multconj(Bk,Zetamc)+Zetamc*Ck);
      h[k]=h[k]*ninv+Bk+Ck;
    }
    if(even) h[c]=h[c]*ninv+multconj(B[c].real(),Zetak)+Zetak*C[c].real();
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
