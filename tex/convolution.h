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
  // Pass a temporary work array f to save memory.
  convolution(unsigned int n, unsigned int m, Complex *f=NULL) :
    n(n), m(m) {
    if(f) {
      F=NULL;
      G=NULL;
      rc=new rcfft1d(n,f);
      cr=new crfft1d(n,f);
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

class cconvolution {
protected:
  unsigned int n;
  unsigned int m;
  fft1d *Backwards;
  fft1d *Forwards;
  double c,s;
public:  
  cconvolution(unsigned int n, unsigned int m, Complex *f) :
    n(n), m(m) {
    Backwards=new fft1d(n,1,f);
    Forwards=new fft1d(n,-1,f);
  }
  
  cconvolution(unsigned int m, Complex *f) : m(m) {
    n=2*m;
    
    double arg=2.0*M_PI/n;
    c=cos(arg);
    s=sin(arg);

    Forwards=new fft1d(m,1,f);
    Backwards=new fft1d(m,-1,f);
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
    for(unsigned int i=m; i < n; i++) f[i]=0.0;
    Backwards->fft(f);
  
    for(unsigned int i=m; i < n; i++) g[i]=0.0;
    Backwards->fft(g);
      
    double ninv=1.0/n;
    for(unsigned int i=0; i < n; ++i)
      f[i] *= g[i]*ninv;
	
    Forwards->fft(f);
  }
  
  // Note: input arrays f and g are destroyed.
  // u is a temporary work array of size n.
  void unpadded(Complex *f, Complex *g, Complex *u) {
    double re=1.0;
    double im=0.0;
    for(unsigned int k=0; k < m; ++k) {
      Complex *P=u+k;
      Complex *Q=P+m;
      Complex fk=*(f+k);
      Complex gk=*(g+k);
      P->re=re*fk.re-im*fk.im;
      P->im=im*fk.re+re*fk.im;
      Q->re=re*gk.re-im*gk.im;
      Q->im=im*gk.re+re*gk.im;
      double temp=re*c+im*s; 
      im=-re*s+im*c;
      re=temp;
    }  
    
    Backwards->fft(f);
    Backwards->fft(u);
    Backwards->fft(g);
    Backwards->fft(u+m);
    
    for(unsigned int k=0; k < m; ++k) {
      Complex *p=f+k;
      Complex fk=*p;
      Complex gk=*(g+k);
      p->re=fk.re*gk.re-fk.im*gk.im;
      p->im=fk.re*gk.im+fk.im*gk.re;
    }
    
    Forwards->fft(f);
    
    for(unsigned int k=0; k < m; ++k) {
      Complex *p=u+k;
      Complex fk=*p;
      Complex gk=*(p+m);
      p->re=fk.re*gk.re-fk.im*gk.im;
      p->im=fk.re*gk.im+fk.im*gk.re;
    }
    
    Forwards->fft(u);
    
    double ninv=1.0/n;
    re=ninv;
    im=0.0;
    for(unsigned int k=0; k < m; ++k) {
      Complex *p=f+k;
      Complex fk=*p;
      Complex fkm=*(u+k);
      p->re=ninv*fk.re+re*fkm.re-im*fkm.im;
      p->im=ninv*fk.im+im*fkm.re+re*fkm.im;
      double temp=re*c-im*s;
      im=re*s+im*c;
      re=temp;
    }
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
      H[i]=sum;
    }
  }	

};

// Untested!
class mcconvolution {
protected:
  unsigned int n;
  unsigned int m;
  unsigned int M;
  unsigned int stride;
  mfft1d *Backwards;
  mfft1d *Forwards;
  double c,s;
public:  
  mcconvolution(unsigned int m, unsigned int M, unsigned int stride,
                Complex *f) : m(m), M(M), stride(stride) {
    n=2*m;
    
    double arg=2.0*M_PI/n;
    c=cos(arg);
    s=sin(arg);

    Backwards=new mfft1d(m,-1,M,stride,1,f);
    Forwards=new mfft1d(m,1,M,stride,1,f);
  }
  
  // Note: input arrays f and g are destroyed.
  // u is a temporary work array of size n*M.
  void unpadded(Complex *f, Complex *g, Complex *u) {
    double re=1.0;
    double im=0.0;
    unsigned int stop=m*stride;
    unsigned int mM=m*M;
    
    for(unsigned int k=0; k < stop; k += stride) {
      Complex *uk=u+k;
      Complex *Fk=f+k;
      Complex *Gk=g+k;
      for(unsigned int i=0; i < M; ++i) {
        Complex *P=uk+i;
        Complex *Q=P+mM;
        Complex fk=*(Fk+i);
        Complex gk=*(Gk+i);
        P->re=re*fk.re-im*fk.im;
        P->im=im*fk.re+re*fk.im;
        Q->re=re*gk.re-im*gk.im;
        Q->im=im*gk.re+re*gk.im;
      }  
      double temp=re*c+im*s; 
      im=-re*s+im*c;
      re=temp;
    }
    
    Forwards->fft(f);
    Forwards->fft(u);
    Forwards->fft(g);
    Forwards->fft(u+mM);
    
    for(unsigned int k=0; k < stop; k += stride) {
      Complex *fk=f+k;
      Complex *Gk=g+k;
      for(unsigned int i=0; i < M; ++i) {
        Complex *p=fk+i;
        Complex fk=*p;
        Complex gk=*(Gk+i);
        p->re=fk.re*gk.re-fk.im*gk.im;
        p->im=fk.re*gk.im+fk.im*gk.re;
      }
    }
    
    Backwards->fft(f);
    
    for(unsigned int k=0; k < stop; k += stride) {
      Complex *uk=u+k;
      for(unsigned int i=0; i < M; ++i) {
        Complex *p=uk+i;
        Complex fk=*p;
        Complex gk=*(p+mM);
        p->re=fk.re*gk.re-fk.im*gk.im;
        p->im=fk.re*gk.im+fk.im*gk.re;
      }
    }
    
    
    Backwards->fft(u);
    
    double ninv=1.0/n;
    re=ninv;
    im=0.0;
    for(unsigned int k=0; k < stop; k += stride) {
      Complex *uk=u+k;
      Complex *fk=f+k;
      for(unsigned int i=0; i < M; ++i) {
        Complex *p=fk+i;
        Complex fk=*p;
        Complex fkm=*(uk+i);
        p->re=ninv*fk.re+re*fkm.re-im*fkm.im;
        p->im=ninv*fk.im+im*fkm.re+re*fkm.im;
      }
      double temp=re*c-im*s;
      im=re*s+im*c;
      re=temp;
    }
  }
};

// Compute the virtual m-padded complex Fourier transform of M complex
// vectors, each of length m.
// Before calling fft(), the arrays in and out (which may coincide) must be
// allocated as Complex[M*m].
//
// In-place usage:
//
//   ffthalf Backward(m,M,stride);
//   Backward.fft(in);
//
// Notes:
//   stride is the spacing between the elements of each Complex vector;
//
class ffthalf {
  unsigned int n;
  unsigned int m;
  unsigned int M;
  unsigned int stride;
  unsigned int dist;
  mfft1d *Backwards;
  mfft1d *Forwards;
  Complex zeta;
  double c,s;

public:  
  ffthalf(unsigned int m, unsigned int M,
          unsigned int stride, Complex *f) : m(m), M(M), stride(stride) {
    n=2*m;
    double arg=2.0*M_PI/n;
    c=cos(arg);
    s=sin(arg);
    
    // TODO: Standardize signs
    Backwards=new mfft1d(m,-1,M,stride,1,f);
    Forwards=new mfft1d(m,1,M,stride,1,f);
  }
  
  void backwards(Complex *f, Complex *u) {
    double re=1.0;
    double im=0.0;
    unsigned int stop=m*stride;
    for(unsigned int k=0; k < stop; k += stride) {
      Complex *fk=f+k;
      Complex *uk=u+k;
      for(unsigned int i=0; i < M; ++i) {
        Complex *P=uk+i;
        Complex *p=fk+i;
        Complex fk=*p;
        P->re=re*fk.re-im*fk.im;
        P->im=im*fk.re+re*fk.im;
      }
      double temp=re*c+im*s; 
      im=-re*s+im*c;
      re=temp;
    }
    
    Backwards->fft(f);
    Backwards->fft(u);
  }
  
  void forwards(Complex *f, Complex *u) {
    Forwards->fft(f);
    Forwards->fft(u);

    double ninv=1.0/n;
    double re=ninv;
    double im=0.0;
    unsigned int stop=m*stride;
    for(unsigned int k=0; k < stop; k += stride) {
      Complex *uk=u+k;
      Complex *fk=f+k;
      for(unsigned int i=0; i < M; ++i) {
        Complex *p=fk+i;
        Complex fk=*p;
        Complex fkm=*(uk+i);
        p->re=ninv*fk.re+re*fkm.re-im*fkm.im;
        p->im=ninv*fk.im+im*fkm.re+re*fkm.im;
      }
      double temp=re*c-im*s;
      im=re*s+im*c;
      re=temp;
    }
  }
};
  
class cconvolution2 {
protected:
  unsigned int n;
  unsigned int m;
  mfft1d *xBackwards;
  mfft1d *yBackwards;
  mfft1d *xForwards;
  mfft1d *yForwards;
  fft2d *Backwards;
  fft2d *Forwards;
  cconvolution *C;
  ffthalf *fftpad;
  Complex *u,*v;
  Complex *work;
//  mcconvolution *MC;
//  Complex *mwork;
public:  
  cconvolution2(unsigned int n, unsigned int m, Complex *f) :
    n(n), m(m) {
    xBackwards=new mfft1d(n,1,m,n,1,f);
    yBackwards=new mfft1d(n,1,n,1,n,f);
    xForwards=new mfft1d(n,-1,m,n,1,f);
    yForwards=new mfft1d(n,-1,n,1,n,f);
    Backwards=new fft2d(n,n,1,f);
    Forwards=new fft2d(n,n,-1,f);
  }
  
  cconvolution2(unsigned int m, Complex *f) : m(m) {
    n=2*m;
    u=FFTWComplex(m*m);
    v=FFTWComplex(m*m);
    work=FFTWComplex(n);
    fftpad=new ffthalf(m,m,m,u);
    C=new cconvolution(m,work);
//    mwork=FFTWComplex(2*m*m);
//    MC=new mcconvolution(m,m,m,mwork);
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

  void pad(Complex *f) {
    for(unsigned int i=0; i < m;) {
      unsigned int j=n*i+m;
      ++i;
      unsigned int nip=n*i;
      for(; j < nip; ++j)
        f[j]=0.0;
    }
    
    for(unsigned int i=m; i < n;) {
      unsigned int j=n*i;
      ++i;
      unsigned int nip=n*i;
      for(; j < nip; ++j)
        f[j]=0.0;
    }
  }
  
  void fft(Complex *f, Complex *g) {
    pad(f);
//    Backwards->fft(f);
    xBackwards->fft(f);
    yBackwards->fft(f);
  
    pad(g);
    Backwards->fft(g);
    
    unsigned int n2=n*n;
    double ninv=1.0/n2;
    for(unsigned int i=0; i < n2; ++i)
        f[i] *= g[i]*ninv;
	
//    Forwards->fft(f);
    yForwards->fft(f);
    xForwards->fft(f);
  }
  
  // Note: input arrays f and g are destroyed.
  void unpadded(Complex *f, Complex *g) {
    fftpad->backwards(f,u);
    fftpad->backwards(g,v);

//    MC->unpadded(f,g,mwork);
//    MC->unpadded(u,v,mwork);
    unsigned int m2=m*m;
    for(unsigned int i=0; i < m2; i += m)
      C->unpadded(f+i,g+i,work);
    for(unsigned int i=0; i < m2; i += m)
      C->unpadded(u+i,v+i,work);
    
    fftpad->forwards(f,u);
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
    for(unsigned int i=0; i < m; ++i) {
      for(unsigned int j=0; j < m; ++j) {
        Complex sum=0.0;
        for(unsigned int k=0; k <= i; ++k)
          for(unsigned int p=0; p <= j; ++p)
            sum += F[k*m+p]*G[(i-k)*m+j-p];
        H[i*m+j]=sum;
      }
    }
  }	

};

#endif
