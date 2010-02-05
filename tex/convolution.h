#include "fftw++.h"

#ifndef __convolution_h__
#define __convolution_h__ 1

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

static const double sqrt3=sqrt(3.0);
static const double hsqrt3=0.5*sqrt3;

class convolution {
protected:
  unsigned int n;
  unsigned int m;
  unsigned int c;
  rcfft1d *rc, *rco;
  crfft1d *cr, *cro;
  double *F,*G;
  double Cos;
  double Sin;
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
  
  convolution(unsigned int m) : m(m) {}
  
  // u and v must be each allocated as m/2+1 Complex values.
  convolution(unsigned int m, Complex *u, Complex *v) : m(m) {
    n=3*m;
    c=m/2;
    double arg=2.0*M_PI/n;
    Cos=cos(arg);
    Sin=sin(arg);

    rc=new rcfft1d(m,u);
    cr=new crfft1d(m,u);

    double *U=(double *) u;
    rco=new rcfft1d(m,U,v);
    cro=new crfft1d(m,v,U);

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
  // u and v are temporary work arrays each of size m/2+1.
  void unpadded(Complex *f, Complex *g, Complex *u, Complex *v) {
    double f0=f[0].re;
    double g0=g[0].re;

    bool even=m % 2 == 0;
    if(!even) _exit(1);
    
    // Arrange input data
    u[0]=f0;
    v[0]=g0;
    Complex fc=f[c];
    unsigned int m1=m-1;
    Complex fmk=f[m1];
    double fmkre=fmk.re;
    double fmkim=fmk.im;
    f[m1]=f0;
    Complex gc=g[c];
    Complex gmk=g[m1];
    double gmkre=gmk.re;
    double gmkim=gmk.im;
    g[m1]=g0;
    
    double Re=Cos;
    double Im=-Sin;
    
    for(unsigned int k=1; k < c; ++k) {
      Complex *p=f+k;
      double re=-0.5*fmkre+p->re;
      double im=hsqrt3*fmkre;
      double Are=Re*re-Im*im;
      double Aim=Re*im+Im*re;
      re=-0.5*fmkim-p->im;
      im=hsqrt3*fmkim;
      p->re += fmkre;
      p->im -= fmkim;
      double Bre=-Re*im-Im*re;
      double Bim=Re*re-Im*im;
      p=u+k;
      p->re=Are-Bre;
      p->im=Aim-Bim;
      p=f+m1-k;
      fmkre=p->re;
      fmkim=p->im;
      p->re=Are+Bre;
      p->im=Aim+Bim;

      p=g+k;
      re=-0.5*gmkre+p->re;
      im=hsqrt3*gmkre;
      Are=Re*re-Im*im;
      Aim=Re*im+Im*re;
      re=-0.5*gmkim-p->im;
      im=hsqrt3*gmkim;
      p->re += gmkre;
      p->im -= gmkim;
      Bre=-Re*im-Im*re;
      Bim=Re*re-Im*im;
      p=v+k;
      p->re=Are-Bre;
      p->im=Aim-Bim;
      p=g+m1-k;
      gmkre=p->re;
      gmkim=p->im;
      p->re=Are+Bre;
      p->im=Aim+Bim;
      double temp=Re*Cos+Im*Sin; 
      Im=-Re*Sin+Im*Cos;
      Re=temp;
    }
  
    double A=fc.re;
    double B=sqrt3*fc.im;
    fc=f[c];
    f[c]=2.0*A;
    u[c]=A+B;
    A -= B;
    
    double C=gc.re;
    B=sqrt3*gc.im;
    gc=g[c];
    g[c]=2.0*C;
    v[c]=C+B;
    C -= B;
    
    // FFTs, convolutions, and inverse FFTs
    // seven of nine transforms are out-of-place
    // r=-1:
    cr->fft(u);
    cr->fft(v);
    double *U=(double *) u;
    double *V=(double *) v;
    for(int i=0; i < m; ++i)
      V[i] *= U[i];
    rco->fft(v,U); // v is now free

    // r=0:
    V=(double *) v;
    cro->fft(f,V);
    double *F=(double *) f;
    cro->fft(g,F);
    for(int i=0; i < m; ++i)
      V[i] *= F[i];
    rco->fft(V,f);
    unsigned int cm1=c-1;
    Complex overlap0=f[cm1];
    double overlap1=f[c].re;

    // r=1:
    f[cm1]=A;
    f[c]=fc;
    Complex *f1=f+cm1;
    Complex *g1=g+cm1;
    G=(double *) g;
    g[cm1]=C;
    g[c]=gc;
    V=(double *) v;
    cro->fft(g1,V);
    cro->fft(f1,G);
    for(int i=0; i < m; ++i)
      G[i] *= V[i];
    rco->fft(G,v);

    unsigned int stop=m-c-1;
    double ninv=1.0/n;
    f[0]=(f[0].re+v[0].re+u[0].re)*ninv;
    Re=Cos*ninv;
    Im=Sin*ninv;
    Complex *fm=f+m;
    for(unsigned k=1; k < stop; ++k) {
      Complex *p=f+k;
      Complex *q=v+k;
      Complex *r=u+k;
      Complex *s=fm-k;
      double f0re=p->re*ninv;
      double f0im=p->im*ninv;
      double f1re=Re*q->re+Im*q->im;
      double f2re=Re*r->re-Im*r->im;
      double sre=f1re+f2re;
      double f1im=Re*q->im-Im*q->re;
      double f2im=Re*r->im+Im*r->re;
      double sim=f1im+f2im;
      p->re=f0re+sre;
      p->im=f0im+sim;
      s->re=f0re-0.5*sre-hsqrt3*(f1im-f2im);
      s->im=-f0im+0.5*sim-hsqrt3*(f1re-f2re);
      double temp=Re*Cos-Im*Sin; 
      Im=Re*Sin+Im*Cos;
      Re=temp;
    }
    
    Complex Zetak=Complex(Re,Im);
    Complex f0k=overlap0*ninv;
    Complex f1k=conj(Zetak)*v[stop];
    Complex f2k=Zetak*u[cm1];
    f[cm1]=f0k+f1k+f2k;
    static const Complex zeta3(-0.5,hsqrt3);
    f[c+1]=conj(f0k+zeta3*f1k)+zeta3*conj(f2k);

    if(even) f[c]=(overlap1-v[c].re*zeta3-u[c].re*conj(zeta3))*ninv;
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


// calculates a the convolution of two complex, non-Hermitian vectors
class cconvolution {
 protected:
  unsigned int n;
  unsigned int m;
  fft1d *Backwards, *Backwardso;
  fft1d *Forwards, *Forwardso;
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

    Backwards=new fft1d(m,-1,f);
    Forwards=new fft1d(m,1,f);
    
    Complex *G=(Complex *)FFTWComplex(m);
    Backwardso=new fft1d(m,-1,f,G);
    Forwardso=new fft1d(m,1,G,f);
    FFTWdelete(G);
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
  // u and v are temporary work arrays each of size m.
  // 1/2 padding
  void unpadded(Complex *f, Complex *g, Complex *u, Complex *v) {
    double re=1.0;
    double im=0.0;
    for(unsigned int k=0; k < m; ++k) {
      Complex *P=u+k;
      Complex *Q=v+k;
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
    
    Backwards->fft(u);
    Backwards->fft(v);
    for(unsigned int k=0; k < m; ++k) {
      Complex *p=v+k;
      Complex vk=*p;
      Complex uk=*(u+k);
      p->re=uk.re*vk.re-uk.im*vk.im;
      p->im=uk.re*vk.im+uk.im*vk.re;
    }
    Forwardso->fft(v,u);

    Backwardso->fft(f,v);
    Backwardso->fft(g,f);
    for(unsigned int k=0; k < m; ++k) {
      Complex *p=v+k;
      Complex vk=*p;
      Complex fk=*(f+k);
      p->re=vk.re*fk.re-vk.im*fk.im;
      p->im=vk.re*fk.im+vk.im*fk.re;
    }
    Forwardso->fft(v,f);
    
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

class mcconvolution {
protected:
  unsigned int n;
  unsigned int m;
  unsigned int M;
  unsigned int stride;
  unsigned int dist;
  mfft1d *Backwards;
  mfft1d *Forwards;
  double c,s;
public:  
  mcconvolution(unsigned int m, unsigned int M, unsigned int stride,
                unsigned int dist, Complex *f) : m(m), M(M), stride(stride),
                                                 dist(dist) {
    n=2*m;
    
    double arg=2.0*M_PI/n;
    c=cos(arg);
    s=sin(arg);

    Backwards=new mfft1d(m,-1,M,stride,dist,f);
    Forwards=new mfft1d(m,1,M,stride,dist,f);
  }
  
  // Special case optimized for stride=1.
  void unpadded1(Complex *f, Complex *g, Complex *u) {
    unsigned int istop=M*dist;
    unsigned int mM=m*M;
    for(unsigned int i=0; i < istop; i += dist) {
      double re=1.0;
      double im=0.0;
      Complex *Gi=g+i;
      Complex *Fi=f+i;
      Complex *ui=u+i;
      for(unsigned int k=0; k < m; ++k) {
        Complex *P=ui+k;
        Complex *Q=P+mM;
        Complex fk=*(Fi+k);
        Complex gk=*(Gi+k);
        P->re=re*fk.re-im*fk.im;
        P->im=im*fk.re+re*fk.im;
        Q->re=re*gk.re-im*gk.im;
        Q->im=im*gk.re+re*gk.im;
        double temp=re*c+im*s; 
        im=-re*s+im*c;
        re=temp;
      }  
    }
    
    Forwards->fft(f);
    Forwards->fft(u);
    Forwards->fft(g);
    Forwards->fft(u+mM);
    
    for(unsigned int i=0; i < istop; i += dist) {
      Complex *fi=f+i;
      Complex *Gi=g+i;
      for(unsigned int k=0; k < m; ++k) {
        Complex *p=fi+k;
        Complex fk=*p;
        Complex gk=*(Gi+k);
        p->re=fk.re*gk.re-fk.im*gk.im;
        p->im=fk.re*gk.im+fk.im*gk.re;
      }
    }
    
    Backwards->fft(f);
    
    for(unsigned int i=0; i < istop; i += dist) {
      Complex *ui=u+i;
      for(unsigned int k=0; k < m; ++k) {
        Complex *p=ui+k;
        Complex fk=*p;
        Complex gk=*(p+mM);
        p->re=fk.re*gk.re-fk.im*gk.im;
        p->im=fk.re*gk.im+fk.im*gk.re;
      }
    }
    
    Backwards->fft(u);
    
    double ninv=1.0/n;
    for(unsigned int i=0; i < istop; i += dist) {
      Complex *ui=u+i;
      Complex *fi=f+i;
      double re=ninv;
      double im=0.0;
      for(unsigned int k=0; k < m; ++k) {
        Complex *p=fi+k;
        Complex fk=*p;
        Complex fkm=*(ui+k);
        p->re=ninv*fk.re+re*fkm.re-im*fkm.im;
        p->im=ninv*fk.im+im*fkm.re+re*fkm.im;
        double temp=re*c-im*s;
        im=re*s+im*c;
        re=temp;
      }
    }
  }

  // Note: input arrays f and g are destroyed.
  // u is a temporary work array of size n*M.
  void unpadded(Complex *f, Complex *g, Complex *u) {
    if(stride == 1) {
      unpadded1(f,g,u);
      return;
    }
    double re=1.0;
    double im=0.0;
    unsigned int kstop=m*stride;
    unsigned int istop=M*dist;
    unsigned int mM=m*M;
    for(unsigned int k=0; k < kstop; k += stride) {
      Complex *Gk=g+k;
      Complex *Fk=f+k;
      Complex *uk=u+k;
      for(unsigned int i=0; i < istop; i += dist) {
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
    
    for(unsigned int k=0; k < kstop; k += stride) {
      Complex *fk=f+k;
      Complex *Gk=g+k;
      for(unsigned int i=0; i < istop; i += dist) {
        Complex *p=fk+i;
        Complex fk=*p;
        Complex gk=*(Gk+i);
        p->re=fk.re*gk.re-fk.im*gk.im;
        p->im=fk.re*gk.im+fk.im*gk.re;
      }
    }
    
    Backwards->fft(f);
    
    for(unsigned int k=0; k < kstop; k += stride) {
      Complex *uk=u+k;
      for(unsigned int i=0; i < istop; i += dist) {
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
    for(unsigned int k=0; k < kstop; k += stride) {
      Complex *uk=u+k;
      Complex *fk=f+k;
      for(unsigned int i=0; i < istop; i += dist) {
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

// Compute the scrambled virtual m-padded complex Fourier transform of M complex
// vectors, each of length m.
// Before calling fft(), the arrays in and out (which may coincide), along
// with the array u, must be allocated as Complex[M*m].
//
//   ffthalf fftpad(m,M,stride);
//   fftpad.backwards(in,u);
//   fftpad.forwards(in,u);
//
// Notes:
//   stride is the spacing between the elements of each Complex vector.
//
class ffthalf {
  unsigned int n;
  unsigned int m;
  unsigned int M;
  unsigned int stride;
  unsigned int dist;
  mfft1d *Backwards;
  mfft1d *Forwards;
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
  
// Compute the scrambled virtual m-padded complex Fourier transform of M complex
// vectors, each of length 2m-1.
// Before calling fft(), the arrays in and out (which may coincide)
// must be allocated as Complex[M*(2m-1)].
// The array u must be allocated as Complex[M*(m+1)].
// The array work must be allocated as Complex[M].
//
//   ffttwothirds fftpad(m,M,stride);
//   fftpad.backwards(in,u);
//   fftpad.forwards(in,u);
//
// Notes:
//   stride is the spacing between the elements of each Complex vector.
//
class ffttwothirds {
  unsigned int n;
  unsigned int m;
  unsigned int M;
  unsigned int stride;
  unsigned int dist;
  mfft1d *Backwards;
  mfft1d *Forwards;
  double Cos,Sin;
  Complex zeta;
public:  
  ffttwothirds(unsigned int m, unsigned int M, unsigned int stride,
               Complex *f) : m(m), M(M), stride(stride) {
    n=3*m;
    double arg=2.0*M_PI/n;
    Cos=cos(arg);
    Sin=sin(arg);
    zeta=Complex(cos(arg),sin(arg));
    
    Backwards=new mfft1d(m,1,M,stride,1,f);
    Forwards=new mfft1d(m,-1,M,stride,1,f);
  }
  
  void backwards(Complex *f, Complex *u) {
    double Re=Cos;
    double Im=Sin;
    Complex fk0=f[0];
    double fkre=fk0.re;
    double fkim=fk0.im;
    unsigned int m1=m-1;
    unsigned int m1stride=m1*stride;
    Complex *fm1stride=f+m1stride;
    for(unsigned int i=0; i < M; ++i) {
      u[i]=f[i]=fm1stride[i];
    }
    
    unsigned int stop=m1stride;
    for(unsigned int k=stride; k < stop; k += stride) {
      Complex *uk=u+k;
      Complex *fk=f+k;
      unsigned int mk=m1+k;
      Complex *fmk=f+mk;
      for(unsigned int i=0; i < M; ++i) {
        Complex *p=fmk+i;
        double fmkre=p->re;
        double re=-0.5*fkre+fmkre;
        double im=-hsqrt3*fkre;
        double Are=Re*re-Im*im;
        double Aim=Re*im+Im*re;
        re=p->im-0.5*fkim;
        im=-hsqrt3*fkim;
        double Bre=-Re*im-Im*re;
        double Bim=Re*re-Im*im;
        re=fkre+fmkre;
        im=fkim+p->im;
        p->re=Are+Bre;
        p->im=Aim+Bim;
        p=fk+i;
        fkre=p->re;
        fkim=p->im;
        p->re=re;
        p->im=im;
        p=uk+i;
        p->re=Are-Bre;
        p->im=Bim-Aim;
      }
      double temp=Re*Cos-Im*Sin;
      Im=Re*Sin+Im*Cos;
      Re=temp;
    }
      
    Complex *uk=u+m1stride;
    Complex *fk=fm1stride;
    unsigned int mk=m1+m1stride;
    Complex *p=f+mk;
    for(unsigned int i=0; i < M; ++i) {
      Complex fmk=p[i];
      Complex Zetak=Complex(Re,Im);
      static Complex zeta3c(-0.5,-hsqrt3);
      Complex A=Zetak*(fmk.re+zeta3c*fkre);
      Complex B=Complex(0,1)*Zetak*(fmk.im+zeta3c*fkim);
      f[m1stride+i]=Complex(fkre,fkim)+fmk;
      f[mk+i]=A+B;
      u[m1stride+i]=conj(A-B);
    }
      
    Backwards->fft(f);
    for(unsigned int i=0; i < M; ++i) {
      u[stride*m+i]=f[stride*(m-1)+i]; // Store extra value here.
      f[stride*(m-1)+i]=u[i];
    }
    
    Backwards->fft(f+(m-1)*stride);
    Backwards->fft(u);
  }
  
  void forwards(Complex *f, Complex *u) {
    unsigned int m1=m-1;
    Complex *fm1=f+m1;
    Forwards->fft(fm1);
    Complex overlap0=*fm1;
    *fm1=u[m];
    Forwards->fft(f);
    Forwards->fft(u);

    double ninv=1.0/n;
    overlap0=(f[0]+overlap0+u[0])*ninv;
    double Re=Cos*ninv;
    double Im=Sin*ninv;
    for(int k=1; k < m; ++k) {
      Complex *p=f+k;
      Complex *q=fm1+k;
      Complex *r=u+k;
      double f0re=p->re*ninv;
      double f0im=p->im*ninv;
      double f1re=Re*q->re+Im*q->im;
      double f1im=Re*q->im-Im*q->re;
      double f2re=Re*r->re-Im*r->im;
      double f2im=Re*r->im+Im*r->re;
      double sre=f1re+f2re;
      double sim=f1im+f2im;
      --p;
      p->re=f0re-0.5*sre-hsqrt3*(f1im-f2im);
      p->im=f0im-0.5*sim+hsqrt3*(f1re-f2re);
      q->re=f0re+sre;
      q->im=f0im+sim;
      double temp=Re*Cos-Im*Sin;
      Im=Re*Sin+Im*Cos;
      Re=temp;
    }
    f[m1]=overlap0;
  }
};
  
class cconvolution2 {
protected:
  unsigned int n;
  unsigned int m;
  bool prune;
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
public:  
  // Set prune=true to skip Fourier transforming zero rows.
  cconvolution2(unsigned int n, unsigned int m, Complex *f, bool prune=false) :
    n(n), m(m), prune(prune) {
    if(prune) {
      xBackwards=new mfft1d(n,1,m,n,1,f);
      yBackwards=new mfft1d(n,1,n,1,n,f);
      xForwards=new mfft1d(n,-1,m,n,1,f);
      yForwards=new mfft1d(n,-1,n,1,n,f);
    } else {
      Backwards=new fft2d(n,n,1,f);
      Forwards=new fft2d(n,n,-1,f);
    }
  }
  
  cconvolution2(unsigned int m, Complex *f) : m(m) {
    n=2*m;
    u=FFTWComplex(m*m);
    v=FFTWComplex(m*m);
    work=FFTWComplex(n);
    fftpad=new ffthalf(m,m,m,u);
    C=new cconvolution(m,work);
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
    if(prune) {
      xBackwards->fft(f);
      yBackwards->fft(f);
    } else
      Backwards->fft(f);
  
    pad(g);
    if(prune) {
      xBackwards->fft(g);
      yBackwards->fft(g);
    } else
      Backwards->fft(g);
    
    unsigned int n2=n*n;
    double ninv=1.0/n2;
    for(unsigned int i=0; i < n2; ++i)
      f[i] *= g[i]*ninv;
	
    if(prune) {
      yForwards->fft(f);
      xForwards->fft(f);
    } else {
      Forwards->fft(f);
    }
  }
  
  // Note: input arrays f and g are destroyed.
  void unpadded(Complex *f, Complex *g) {
    fftpad->backwards(f,u);
    fftpad->backwards(g,v);

    unsigned int m2=m*m;
    Complex *work2=work+m;
    for(unsigned int i=0; i < m2; i += m)
      C->unpadded(f+i,g+i,work,work2);
    for(unsigned int i=0; i < m2; i += m)
      C->unpadded(u+i,v+i,work,work2);
    
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

#if 0
class convolution2 {
protected:
  unsigned int n;
  unsigned int m;
  bool prune;
  mfft1d *xBackwards;
  mfft1d *yBackwards;
  mfft1d *xForwards;
  mfft1d *yForwards;
  fft2d *Backwards;
  fft2d *Forwards;
  convolution *C;
  ffttwothirds *fftpad;
  Complex *u,*v;
  Complex *work;
public:  
  // Set prune=true to skip Fourier transforming zero rows.
  convolution2(unsigned int n, unsigned int m, Complex *f, bool prune=false) :
    n(n), m(m), prune(prune) {
    if(prune) {
      xBackwards=new mfft1d(n,1,m,n,1,f);
      yBackwards=new mfft1d(n,1,n,1,n,f);
      xForwards=new mfft1d(n,-1,m,n,1,f);
      yForwards=new mfft1d(n,-1,n,1,n,f);
    } else {
      Backwards=new fft2d(n,n,1,f);
      Forwards=new fft2d(n,n,-1,f);
    }
  }
  
  convolution2(unsigned int m, Complex *f) : m(m) {
    n=3*m;
    u=FFTWComplex(m*m);
    v=FFTWComplex(m*m);
    work=FFTWComplex(n);
    fftpad=new ffthalf(m,m,m,u);
    C=new cconvolution(m,work);
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
  }
  
  void fft(Complex *f, Complex *g) {
    pad(f);
    if(prune) {
      xBackwards->fft(f);
      yBackwards->fft(f);
    } else
      Backwards->fft(f);
  
    pad(g);
    if(prune) {
      xBackwards->fft(g);
      yBackwards->fft(g);
    } else
      Backwards->fft(g);
    
    unsigned int n2=n*n;
    double ninv=1.0/n2;
    for(unsigned int i=0; i < n2; ++i)
      f[i] *= g[i]*ninv;
	
    if(prune) {
      yForwards->fft(f);
      xForwards->fft(f);
    } else {
      Forwards->fft(f);
    }
  }
  
  // Note: input arrays f and g are destroyed.
  void unpadded(Complex *f, Complex *g) {
    fftpad->backwards(f,u);
    fftpad->backwards(g,v);

    unsigned int m2=m*m;
    Complex *work2=work+m;
    for(unsigned int i=0; i < m2; i += m)
      C->unpadded(f+i,g+i,work,work2);
    for(unsigned int i=0; i < m2; i += m)
      C->unpadded(u+i,v+i,work,work2);
    
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

#endif
