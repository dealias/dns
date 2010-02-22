// TODO: try using simd instructions everywhere for the real term-by-term
// multiply.

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
extern const Complex one;

// In-place explicitly dealiased complex convolution.
class ExplicitConvolution {
protected:
  unsigned int n,m;
  fft1d *Backwards,*Forwards;
  double Cos,Sin;
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
  
  // The distinct) input arrays f and g are each of size n (contents not
  // preserved). The output is returned in f.
  void convolve(Complex *f, Complex *g) {
    for(unsigned int i=m; i < n; i++) f[i]=0.0;
    Backwards->fft(f);
  
    for(unsigned int i=m; i < n; i++) g[i]=0.0;
    Backwards->fft(g);
      
    double ninv=1.0/n;
    for(unsigned int i=0; i < n; ++i)
      f[i] *= g[i]*ninv;
	
    Forwards->fft(f);
  }
};

// In-place implicitly dealiased complex convolution.
class ImplicitConvolution {
protected:
  unsigned int n,m;
  fft1d *Backwards,*Forwards;
  fft1d *Backwardso,*Forwardso;
  double Cos,Sin;
public:  
  
  // u and v are distinct temporary arrays each of size m.
  ImplicitConvolution(unsigned int m, Complex *u, Complex *v) : n(2*m), m(m) {
    double arg=2.0*M_PI/n;
    Cos=cos(arg);
    Sin=sin(arg);

    Backwards=new fft1d(m,1,u);
    Forwards=new fft1d(m,-1,u);
    
    Backwardso=new fft1d(m,1,u,v);
    Forwardso=new fft1d(m,-1,u,v);
  }
  
  ~ImplicitConvolution() {
    delete Forwardso;
    delete Backwardso;
    delete Forwards;
    delete Backwards;
  }
  
  // In-place implicitly dealiased convolution.
  // The input arrays f and g are each of size m (contents not preserved).
  // The output is returned in f.
  // u and v are temporary arrays each of size m.
  void convolve(Complex *f, Complex *g, Complex *u, Complex *v) {
#ifdef __SSE2__      
    const Complex cc(Cos,Cos);
    const Complex ss(-Sin,Sin);
    Vec CC=LOAD(&cc);
    Vec SS=LOAD(&ss);
    Vec Zetak=LOAD(&one);
#else    
    double re=1.0;
    double im=0.0;
#endif
    for(unsigned int k=0; k < m; ++k) {
#ifdef __SSE2__      
      STORE(u+k,ZMULT(Zetak,LOAD(f+k)));
      STORE(v+k,ZMULT(Zetak,LOAD(g+k)));
      Zetak=ZMULT(CC,SS,Zetak);
#else      
      Complex *P=u+k;
      Complex *Q=v+k;
      Complex fk=*(f+k);
      Complex gk=*(g+k);
      P->re=re*fk.re-im*fk.im;
      P->im=im*fk.re+re*fk.im;
      Q->re=re*gk.re-im*gk.im;
      Q->im=im*gk.re+re*gk.im;
      double temp=re*Cos-im*Sin; 
      im=re*Sin+im*Cos;
      re=temp;
#endif      
    }  
    
    // four of six FFTs are out-of-place
    
    Backwards->fft(u);
    Backwards->fft(v);
    for(unsigned int k=0; k < m; ++k) {
#ifdef __SSE2__      
      STORE(v+k,ZMULT(LOAD(u+k),LOAD(v+k)));
#else      
      Complex *p=v+k;
      Complex vk=*p;
      Complex uk=*(u+k);
      p->re=uk.re*vk.re-uk.im*vk.im;
      p->im=uk.re*vk.im+uk.im*vk.re;
#endif      
    }
    Forwardso->fft(v,u);

    Backwardso->fft(f,v);
    Backwardso->fft(g,f);
    for(unsigned int k=0; k < m; ++k) {
#ifdef __SSE2__      
      STORE(v+k,ZMULT(LOAD(v+k),LOAD(f+k)));
#else      
      Complex *p=v+k;
      Complex vk=*p;
      Complex fk=*(f+k);
      p->re=vk.re*fk.re-vk.im*fk.im;
      p->im=vk.re*fk.im+vk.im*fk.re;
#endif      
    }
    Forwardso->fft(v,f);
    
    double ninv=1.0/n;
#ifdef __SSE2__      
    SS=-LOAD(&ss);
    const Complex Ninv(ninv,0.0);
    const Complex Ninv2(ninv,ninv);
    Zetak=LOAD(&Ninv);
    Vec ninv2=LOAD(&Ninv2);
#else    
    re=ninv;
    im=0.0;
#endif    
    for(unsigned int k=0; k < m; ++k) {
#ifdef __SSE2__
      STORE(f+k,ZMULT(Zetak,LOAD(u+k))+ninv2*LOAD(f+k));
      Zetak=ZMULT(CC,SS,Zetak);
#else      
      Complex *p=f+k;
      Complex fk=*p;
      Complex fkm=*(u+k);
      p->re=ninv*fk.re+re*fkm.re-im*fkm.im;
      p->im=ninv*fk.im+im*fkm.re+re*fkm.im;
      double temp=re*Cos+im*Sin;
      im=-re*Sin+im*Cos;
      re=temp;
#endif      
    }
  }
};

// Out-of-place direct complex convolution.
class DirectConvolution {
protected:
  unsigned int m;
public:  
  DirectConvolution(unsigned int m) : m(m) {}
  
  void convolve(Complex *h, Complex *f, Complex *g) {
    for(unsigned int i=0; i < m; i++) {
      Complex sum=0.0;
      for(unsigned int j=0; j <= i; j++) sum += f[j]*g[i-j];
      h[i]=sum;
    }
  }	
};

// In-place explicitly dealiased Hermitian convolution.
class ExplicitHConvolution {
protected:
  unsigned int n;
  unsigned int m;
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
    
  void pad(Complex *f) {
    unsigned int n2=n/2;
    for(unsigned int i=m; i <= n2; i++) f[i]=0.0;
  }
  
// Compute f (*) g, where f and g contain the m non-negative Fourier
// components of real functions. Dealiasing is internally implemented via
// explicit zero-padding to size n >= 3*m.
//
// The (distinct) input arrays f and g must each be allocated to size n/2+1
// (contents not preserved). The output is returned in the first m elements
// of f.
  void convolve(Complex *f, Complex *g) {
    pad(f);
    cr->fft(f);
    
    pad(g);
    cr->fft(g);
	
    double *F=(double *) f;
    double *G=(double *) g;
    double ninv=1.0/n;
    for(unsigned int i=0; i < n; ++i)
      F[i] *= G[i]*ninv;
    rc->fft(f);
  }
};

// In-place implicitly dealiased Hermitian convolution.
class ImplicitHConvolution {
protected:
  unsigned int n;
  unsigned int m;
  unsigned int c;
  rcfft1d *rc, *rco;
  crfft1d *cr, *cro;
  double *F,*G;
  double Cos,Sin;
public:  
  
  // u and v must be each allocated as m/2+1 Complex values.
  ImplicitHConvolution(unsigned int m, Complex *u, Complex *v) : m(m) {
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
  
  ~ImplicitHConvolution() {
    delete cro;
    delete rco;
    delete cr;
    delete rc;
  }
  
  // Note: input arrays f and g are destroyed.
  // f and g are each of size m.
  // u and v are temporary arrays each of size m/2+1.
  void convolve(Complex *f, Complex *g, Complex *u, Complex *v) {
    double f0=f[0].re;
    double g0=g[0].re;

    bool even=m % 2 == 0;
    if(m <= 2 || !even) {
      cout << "Not yet implemented!" << endl;
      _exit(1);
    }
    
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
    
#ifdef __SSE2__      
    const Complex cc(Cos,Cos);
    const Complex ss(-Sin,Sin);
    const Complex zetac(Cos,-Sin);
    Vec CC=LOAD(&cc);
    Vec SS=-LOAD(&ss);
    Vec Zetak=LOAD(&zetac);
    Vec Fmk=LOAD(&fmk);
    Vec Gmk=LOAD(&gmk);
    Vec Mhalf=LOAD(&mhalf);
    Vec HSqrt3=LOAD(&hSqrt3);
#else
    double Re=Cos;
    double Im=-Sin;
#endif
    for(unsigned int k=1; k < c; ++k) {
      Complex *p=f+k;
      Complex *q=g+k;
#ifdef __SSE2__
      Vec A=LOAD(p);
      Vec B=LOAD(q);
      Vec C=Fmk*Mhalf+CONJ(A);
      Vec D=Gmk*Mhalf+CONJ(B);
      STORE(p,A+CONJ(Fmk));
      STORE(q,B+CONJ(Gmk));
      Fmk *= HSqrt3;
      Gmk *= HSqrt3;
      A=ZMULT(Zetak,UNPACKL(C,Fmk));
      B=ZMULTI(Zetak,UNPACKH(C,Fmk));
      C=ZMULT(Zetak,UNPACKL(D,Gmk));
      D=ZMULTI(Zetak,UNPACKH(D,Gmk));
      STORE(u+k,A-B);
      STORE(v+k,C-D);
      p=f+m1-k;
      Fmk=LOAD(p);
      STORE(p,A+B);
      q=g+m1-k;
      Gmk=LOAD(q);
      STORE(q,C+D);
      Zetak=ZMULT(CC,SS,Zetak);
#else
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

      re=-0.5*gmkre+q->re;
      im=hsqrt3*gmkre;
      Are=Re*re-Im*im;
      Aim=Re*im+Im*re;
      re=-0.5*gmkim-q->im;
      im=hsqrt3*gmkim;
      q->re += gmkre;
      q->im -= gmkim;
      Bre=-Re*im-Im*re;
      Bim=Re*re-Im*im;
      q=v+k;
      q->re=Are-Bre;
      q->im=Aim-Bim;
      q=g+m1-k;
      gmkre=q->re;
      gmkim=q->im;
      q->re=Are+Bre;
      q->im=Aim+Bim;
      
      double temp=Re*Cos+Im*Sin; 
      Im=-Re*Sin+Im*Cos;
      Re=temp;
#endif      
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
    
    // seven of nine FFTs are out-of-place
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
    Complex *fm=f+m;
#ifdef __SSE2__      
    Complex Zeta(Cos,Sin);
    const Complex Ninv2(ninv,ninv);
    Vec ninv2=LOAD(&Ninv2);
    Zetak=LOAD(&Zeta)*ninv2;
    SS=LOAD(&ss);
#else
    Re=Cos*ninv;
    Im=Sin*ninv;
#endif    
    for(unsigned k=1; k < stop; ++k) {
      Complex *p=f+k;
      Complex *s=fm-k;
#ifdef __SSE2__      
      Vec F0=LOAD(p)*ninv2;
      Vec F1=ZMULT(CONJ(Zetak),LOAD(v+k));
      Vec F2=ZMULT(Zetak,LOAD(u+k));
      Vec S=F1+F2;
      STORE(p,F0+S);
      STORE(s,CONJ(F0+Mhalf*S)-HSqrt3*FLIP(F1-F2));
      Zetak=ZMULT(CC,SS,Zetak);
#else
      Complex *q=v+k;
      Complex *r=u+k;
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
#endif      
    }
    
#ifdef __SSE2__      
    Complex Zetak0;
    STORE(&Zetak0,Zetak);
#else    
    Complex Zetak0=Complex(Re,Im);
#endif      
    Complex f0k=overlap0*ninv;
    Complex f1k=conj(Zetak0)*v[stop];
    Complex f2k=Zetak0*u[cm1];
    f[cm1]=f0k+f1k+f2k;
    f[c+1]=conj(f0k+zeta3*f1k)+zeta3*conj(f2k);

    if(even) f[c]=(overlap1-v[c].re*zeta3-u[c].re*conj(zeta3))*ninv;
  }
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
  void convolve(Complex *h, Complex *f, Complex *g) {
    for(unsigned int i=0; i < m; i++) {
      Complex sum=0.0;
      for(unsigned int j=0; j <= i; j++) sum += f[j]*g[i-j];
      for(unsigned int j=i+1; j < m; j++) sum += f[j]*conj(g[j-i]);
      for(unsigned int j=1; j < m-i; j++) sum += conj(f[j])*g[i+j];
      h[i]=sum;
    }
  }	
};

// Compute the scrambled virtual m-padded complex Fourier transform of M complex
// vectors, each of length m.
// Before calling fft(), the arrays in and out (which may coincide), along
// with the array u, must be allocated as Complex[M*m].
//
//   fftpad fft(m,M,stride);
//   fft.backwards(in,u);
//   fft.forwards(in,u);
//
// Notes:
//   stride is the spacing between the elements of each Complex vector.
//
class fftpad {
  unsigned int n;
  unsigned int m;
  unsigned int M;
  unsigned int stride;
  unsigned int dist;
  mfft1d *Backwards;
  mfft1d *Forwards;
  double Cos,Sin;

public:  
  fftpad(unsigned int m, unsigned int M,
         unsigned int stride, Complex *f) : m(m), M(M), stride(stride) {
    n=2*m;
    double arg=2.0*M_PI/n;
    Cos=cos(arg);
    Sin=sin(arg);
    
    Backwards=new mfft1d(m,1,M,stride,1,f);
    Forwards=new mfft1d(m,-1,M,stride,1,f);
  }
  
  ~fftpad() {
    delete Forwards;
    delete Backwards;
  }
  
  void backwards(Complex *f, Complex *u) {
#ifdef __SSE2__
    const Complex cc(Cos,Cos);
    const Complex ss(-Sin,Sin);
    const Complex one(1.0,0.0);
    Vec CC=LOAD(&cc);
    Vec SS=LOAD(&ss);
    Vec Zetak=LOAD(&one);
#else
    double re=1.0;
    double im=0.0;
#endif    
    unsigned int stop=m*stride;
    for(unsigned int k=0; k < stop; k += stride) {
      Complex *fk=f+k;
      Complex *uk=u+k;
#ifdef __SSE2__      
      Vec X=UNPACKL(Zetak,Zetak);
      Vec Y=UNPACKH(CONJ(Zetak),Zetak);
      for(unsigned int i=0; i < M; ++i)
        STORE(uk+i,ZMULT(X,Y,LOAD(fk+i)));
      Zetak=ZMULT(CC,SS,Zetak);
#else
      for(unsigned int i=0; i < M; ++i) {
        Complex *p=uk+i;
        Complex fki=*(fk+i);
        p->re=re*fki.re-im*fki.im;
        p->im=im*fki.re+re*fki.im;
      }
      double temp=re*Cos-im*Sin; 
      im=re*Sin+im*Cos;
      re=temp;
#endif      
    }
    
    Backwards->fft(f);
    Backwards->fft(u);
  }
  
  void forwards(Complex *f, Complex *u) {
    Forwards->fft(f);
    Forwards->fft(u);

    double ninv=1.0/n;
#ifdef __SSE2__
    const Complex cc(Cos,Cos);
    const Complex ss(-Sin,Sin);
    const Complex ninv1(ninv,0.0);
    const Complex ninv2(ninv,ninv);
    Vec CC=LOAD(&cc);
    Vec SS=-LOAD(&ss);
    Vec Zetak=LOAD(&ninv1);
    Vec Ninv2=LOAD(&ninv2);
#else
    double re=ninv;
    double im=0.0;
#endif    
    unsigned int stop=m*stride;
    for(unsigned int k=0; k < stop; k += stride) {
      Complex *uk=u+k;
      Complex *fk=f+k;
#ifdef __SSE2__      
      Vec X=UNPACKL(Zetak,Zetak);
      Vec Y=UNPACKH(CONJ(Zetak),Zetak);
      for(unsigned int i=0; i < M; ++i) {
        Complex *p=fk+i;
        STORE(p,LOAD(p)*Ninv2+ZMULT(X,Y,LOAD(uk+i)));
      }
      Zetak=ZMULT(CC,SS,Zetak);
#else        
      for(unsigned int i=0; i < M; ++i) {
        Complex *p=fk+i;
        Complex fki=*p;
        Complex fkm=*(uk+i);
        p->re=ninv*fki.re+re*fkm.re-im*fkm.im;
        p->im=ninv*fki.im+im*fkm.re+re*fkm.im;
      }
      double temp=re*Cos+im*Sin;
      im=-re*Sin+im*Cos;
      re=temp;
#endif     
    }
  }
};
  
// Compute the scrambled virtual m-padded complex Fourier transform of M complex
// vectors, each of length 2m-1 with the origin at index m-1
// (i.e. physical wavenumber k=-m+1 to k=m-1).
// Before calling fft(), the arrays in and out (which may coincide)
// must be allocated as Complex[M*(2m-1)].
// The array u must be allocated as Complex[M*(m+1)].
//
//   fft0pad fft(m,M,stride);
//   fft.backwards(in,u);
//   fft.forwards(in,u);
//
// Notes:
//   stride is the spacing between the elements of each Complex vector.
//
class fft0pad {
  unsigned int n;
  unsigned int m;
  unsigned int M;
  unsigned int stride;
  unsigned int dist;
  mfft1d *Backwards;
  mfft1d *Forwards;
  double Cos,Sin;
public:  
  fft0pad(unsigned int m, unsigned int M, unsigned int stride, Complex *u)
    : m(m), M(M), stride(stride) {
    n=3*m;
    double arg=2.0*M_PI/n;
    Cos=cos(arg);
    Sin=sin(arg);
    
    Backwards=new mfft1d(m,1,M,stride,1,u);
    Forwards=new mfft1d(m,-1,M,stride,1,u);
  }
  
  ~fft0pad() {
    delete Forwards;
    delete Backwards;
  }
  
  void backwards(Complex *f, Complex *u) {
    unsigned int m1=m-1;
    unsigned int m1stride=m1*stride;
    Complex *fm1stride=f+m1stride;
    for(unsigned int i=0; i < M; ++i)
      u[i]=fm1stride[i];
    
    unsigned int mstride=m*stride;
#ifdef __SSE2__      
    const Complex cc(Cos,Cos);
    const Complex ss(-Sin,Sin);
    const Complex zeta(Cos,Sin);
    Vec CC=LOAD(&cc);
    Vec SS=LOAD(&ss);
    Vec Zetak=LOAD(&zeta);
    Vec Mhalf=LOAD(&mhalf);
    Vec Mhsqrt3=LOAD(&mhsqrt3);
#else
    double Re=Cos;
    double Im=Sin;
#endif    
    for(unsigned int k=stride; k < mstride; k += stride) {
      Complex *uk=u+k;
      Complex *fk=f+k;
      Complex *fmk=f+m1stride+k;
#ifdef __SSE2__
      Vec X=UNPACKL(Zetak,Zetak);
      Vec Y=UNPACKH(CONJ(Zetak),Zetak);
      for(unsigned int i=0; i < M; ++i) {
        Complex *p=fmk+i;
        Complex *q=f+i;
        Complex *r=fk+i;
        Vec A=LOAD(p);
        Vec B=LOAD(q);
        Vec Z=B*Mhalf+A;
        STORE(q,LOAD(r));
        STORE(r,B+A);
        B *= Mhsqrt3;
        A=ZMULT(X,Y,UNPACKL(Z,B));
        B=ZMULTI(X,Y,UNPACKH(Z,B));
        STORE(p,A+B);
        STORE(uk+i,CONJ(A-B));
      }
      Zetak=ZMULT(CC,SS,Zetak);
#else        
      for(unsigned int i=0; i < M; ++i) {
        Complex *p=fmk+i;
        Complex *q=f+i;
        double fkre=q->re;
        double fkim=q->im;
        double fmkre=p->re;
        double fmkim=p->im;
        double re=-0.5*fkre+fmkre;
        double im=-hsqrt3*fkre;
        double Are=Re*re-Im*im;
        double Aim=Re*im+Im*re;
        re=fmkim-0.5*fkim;
        im=-hsqrt3*fkim;
        double Bre=-Re*im-Im*re;
        double Bim=Re*re-Im*im;
        p->re=Are+Bre;
        p->im=Aim+Bim;
        p=uk+i;
        p->re=Are-Bre;
        p->im=Bim-Aim;
        p=fk+i;
        q->re=p->re;
        q->im=p->im;
        p->re=fkre+fmkre;
        p->im=fkim+fmkim;
      }
      double temp=Re*Cos-Im*Sin;
      Im=Re*Sin+Im*Cos;
      Re=temp;
#endif      
    }
    
    Backwards->fft(f);
    Complex *umstride=u+mstride;
    for(unsigned int i=0; i < M; ++i) {
      umstride[i]=fm1stride[i]; // Store extra value here.
      fm1stride[i]=u[i];
    }
    
    Backwards->fft(fm1stride);
    Backwards->fft(u);
  }
  
  void forwards(Complex *f, Complex *u) {
    unsigned int m1stride=(m-1)*stride;
    Complex *fm1stride=f+m1stride;
    Forwards->fft(fm1stride);
    unsigned int mstride=m1stride+stride;
    Complex *umstride=u+mstride;
    for(unsigned int i=0; i < M; ++i) {
      Complex temp=umstride[i];
      umstride[i]=fm1stride[i];
      fm1stride[i]=temp;
    }
    
    Forwards->fft(f);
    Forwards->fft(u);

    double ninv=1.0/n;
    for(unsigned int i=0; i < M; ++i)
      umstride[i]=(umstride[i]+f[i]+u[i])*ninv;
#ifdef __SSE2__      
    const Complex cc(Cos,Cos);
    const Complex ss(-Sin,Sin);
    const Complex zeta(Cos,Sin);
    Vec CC=LOAD(&cc);
    Vec SS=LOAD(&ss);
    const Complex Ninv2(ninv,ninv);
    Vec ninv2=LOAD(&Ninv2);
    Vec Zetak=LOAD(&zeta)*ninv2;
    Vec Mhalf=LOAD(&mhalf);
    Vec HSqrt3=LOAD(&hSqrt3);
#else
    double Re=Cos*ninv;
    double Im=Sin*ninv;
#endif    
    for(unsigned int k=stride; k < mstride; k += stride) {
      Complex *fk=f+k;
      Complex *fm1k=fm1stride+k;
      Complex *uk=u+k;
#ifdef __SSE2__      
      Vec X=UNPACKL(Zetak,Zetak);
      Vec Y=UNPACKH(CONJ(Zetak),Zetak);
      for(unsigned int i=0; i < M; ++i) {
        Complex *p=fk+i;
        Complex *q=fm1k+i;
        Vec F0=LOAD(p)*ninv2;
        Vec F1=ZMULT(X,-Y,LOAD(q));
        Vec F2=ZMULT(X,Y,LOAD(uk+i));
        Vec S=F1+F2;
        STORE(p-stride,F0+Mhalf*S+HSqrt3*ZMULTI(F1-F2));
        STORE(q,F0+S);
      }
      Zetak=ZMULT(CC,SS,Zetak);
#else
      for(unsigned int i=0; i < M; ++i) {
        Complex *p=fk+i;
        Complex *q=fm1k+i;
        Complex *r=uk+i;
        double f0re=p->re*ninv;
        double f0im=p->im*ninv;
        double f1re=Re*q->re+Im*q->im;
        double f1im=Re*q->im-Im*q->re;
        double f2re=Re*r->re-Im*r->im;
        double f2im=Re*r->im+Im*r->re;
        double sre=f1re+f2re;
        double sim=f1im+f2im;
        p -= stride;
        p->re=f0re-0.5*sre-hsqrt3*(f1im-f2im);
        p->im=f0im-0.5*sim+hsqrt3*(f1re-f2re);
        q->re=f0re+sre;
        q->im=f0im+sim;
      }
      double temp=Re*Cos-Im*Sin;
      Im=Re*Sin+Im*Cos;
      Re=temp;
#endif      
    }
    for(unsigned int i=0; i < M; ++i)
      fm1stride[i]=umstride[i];
  }
};
  
// In-place explicitly dealiased 2D complex convolution.
class ExplicitConvolution2 {
protected:
  unsigned int n;
  unsigned int m;
  bool prune; // Skip Fourier transforming rows containing all zeroes?
  mfft1d *xBackwards, *xForwards;
  mfft1d *yBackwards, *yForwards;
  fft2d *Backwards, *Forwards;
public:
  ExplicitConvolution2(unsigned int n, unsigned int m,
                       Complex *f, bool prune=false) :
    n(n), m(m), prune(prune) {
    if(prune) {
      xBackwards=new mfft1d(n,1,m,n,1,f);
      yBackwards=new mfft1d(n,1,n,1,n,f);
      yForwards=new mfft1d(n,-1,n,1,n,f);
      xForwards=new mfft1d(n,-1,m,n,1,f);
    } else {
      Backwards=new fft2d(n,n,1,f);
      Forwards=new fft2d(n,n,-1,f);
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
  
  void convolve(Complex *f, Complex *g) {
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
};

// In-place implicitly dealiased 2D complex convolution.
class ImplicitConvolution2 {
protected:
  unsigned int n;
  unsigned int m;
  fftpad *xfftpad;
  ImplicitConvolution *yconvolve;
public:  
  // u1 is a temporary array of size m.
  // u2 is a temporary array of size m*m.
  ImplicitConvolution2(unsigned int m, Complex *u1, Complex *u2) : 
    n(2*m), m(m) {
    xfftpad=new fftpad(m,m,m,u2);
    yconvolve=new ImplicitConvolution(m,u1,u2);
  }
  
  ~ImplicitConvolution2() {
    delete yconvolve;
    delete xfftpad;
  }
  
  // f must be allocated as m*m.
  // u1 and v1 are temporary arrays of size m.
  // u2 and v2 are temporary arrays of size m*m.
  void convolve(Complex *f, Complex *g, Complex *u1, Complex *v1,
                Complex *u2, Complex *v2) {
    xfftpad->backwards(f,u2);
    xfftpad->backwards(g,v2);

    unsigned int m2=m*m;
    for(unsigned int i=0; i < m2; i += m)
      yconvolve->convolve(f+i,g+i,u1,v1);
    for(unsigned int i=0; i < m2; i += m)
      yconvolve->convolve(u2+i,v2+i,u1,v1);
    
    xfftpad->forwards(f,u2);
  }
};

// Out-of-place direct 2D complex convolution.
class DirectConvolution2 {
protected:  
  unsigned int m;
  
public:
  DirectConvolution2(unsigned int m) : m(m) {}
  
  void convolve(Complex *h, Complex *f, Complex *g) {
    for(unsigned int i=0; i < m; ++i) {
      for(unsigned int j=0; j < m; ++j) {
        Complex sum=0.0;
        for(unsigned int k=0; k <= i; ++k)
          for(unsigned int p=0; p <= j; ++p)
            sum += f[k*m+p]*g[(i-k)*m+j-p];
        h[i*m+j]=sum;
      }
    }
  }	
};

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
    if(nx % 2 == 1) prune=true;
    if(prune) {
      xBackwards=new mfft1d(nx,1,my,nyp,1,f);
      yBackwards=new mcrfft1d(ny,nx,1,nyp,f);
      yForwards=new mrcfft1d(ny,nx,1,nyp,f);
      xForwards=new mfft1d(nx,-1,my,nyp,1,f);
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
  
  void pad(Complex *f) {
    unsigned int nyp=ny/2+1;
    unsigned int nx2=nx/2;
    unsigned int end=nx2-mx;
    for(unsigned int i=0; i <= end;) {
      unsigned int j=nyp*i;
      ++i;
      unsigned int nypip=nyp*i;
      for(; j < nypip; ++j)
        f[j]=0.0;
    }
    
    for(unsigned int i=nx2+mx; i < nx;) {
      unsigned int j=nyp*i;
      ++i;
      unsigned int nypip=nyp*i;
      for(; j < nypip; ++j)
        f[j]=0.0;
    }
    for(unsigned int i=0; i < nx;) {
      unsigned int j=nyp*i+my;
      ++i;
      unsigned int nypip=nyp*i;
      for(; j < nypip; ++j)
        f[j]=0.0;
    }
  }
  
  void convolve(Complex *f, Complex *g) {
    pad(f);
    
    if(prune) {
      xBackwards->fft(f);
      fftw::Shift(f,nx,ny,-1);
      yBackwards->fft(f);
    } else
      Backwards->fft(f);
    
    pad(g);
    if(prune) {
      xBackwards->fft(g);
      fftw::Shift(g,nx,ny,-1);
      yBackwards->fft(g);
    } else
      Backwards->fft(g);
    
    double *F=(double *) f;
    double *G=(double *) g;
    
    double ninv=1.0/(nx*ny);
    unsigned int nyp=ny/2+1;
    unsigned int nyp2=2*nyp;

    for(unsigned int i=0; i < nx; ++i) {
      unsigned int nyp2i=nyp2*i;
      unsigned int stop=nyp2i+ny;
      for(unsigned int j=nyp2i; j < stop; ++j)
        F[j] *= G[j]*ninv;
    }
	
    if(prune) {
      yForwards->fft(f);
      fftw::Shift(f,nx,ny,1);
      xForwards->fft(f);
    } else
      Forwards->fft0(f);
  }
};

// In-place implicitly dealiased 2D Hermitian convolution.
class ImplicitHConvolution2 {
protected:
  unsigned int nx,ny;
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

#endif
