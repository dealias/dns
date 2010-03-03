// TODO: auto-symmetrize 2D Hermitian data.
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
extern const Complex one;

static const unsigned int tableSize=1000;

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
    for(unsigned int k=m; k < n; k++) f[k]=0.0;
    Backwards->fft(f);
  
    for(unsigned int k=m; k < n; k++) g[k]=0.0;
    Backwards->fft(g);
      
    double ninv=1.0/n;
#ifdef __SSE2__      
    const Complex ninv2(ninv,ninv);
    Vec Ninv=LOAD(&ninv2);
    for(unsigned int k=0; k < n; ++k)
      STORE(f+k,Ninv*ZMULT(LOAD(f+k),LOAD(g+k)));
#else    
    for(unsigned int k=0; k < n; ++k)
      f[k] *= g[k]*ninv;
#endif    
	
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
#ifdef __SSE2__
  Complex *HalfSec;
#else
  double *HalfSec;
#endif
  Complex *Zeta;
  Complex *Zeta0;
  Complex *ZetaH;
  Complex *ZetaL;
  unsigned int L;
  unsigned int s;
public:  
  
  // u and v are distinct temporary arrays each of size m.
  ImplicitConvolution(unsigned int m, Complex *u, Complex *v) : n(2*m), m(m) {
    double arg=M_PI/m;
    Cos=cos(arg);
    Sin=sin(arg);

    Backwards=new fft1d(m,1,u);
    Forwards=new fft1d(m,-1,u);
    
    Backwardso=new fft1d(m,1,u,v);
    Forwardso=new fft1d(m,-1,u,v);
    
#define SHORTTABLE 0
#if SHORTTABLE    
    s=sqrt(m);
    unsigned int t=(int) ceil(m/s)+1;
    ZetaH=FFTWComplex(t);
    ZetaL=FFTWComplex(s);
    for(unsigned int a=0; a < t; ++a)
      ZetaH[a]=Complex(cos(s*a*arg),sin(s*a*arg));
    for(unsigned int b=0; b < s; ++b)
      ZetaL[b]=Complex(cos(b*arg),sin(b*arg));
#endif
    Zeta=FFTWComplex(m);
    for(unsigned int k=0; k < m; ++k)
      Zeta[k]=Complex(cos(k*arg),sin(k*arg));
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
    Vec Zetak=LOAD(&one);
#else    
    double re=1.0;
    double im=0.0;
#endif
#if SHORTTABLE    
    for(unsigned int a=0, k=0; k < m; ++a) {
      unsigned int stop=min(k+s,m);
      Vec Zeta=LOAD(ZetaH+a);
      Vec X=UNPACKL(Zeta,Zeta);
      Vec Y=UNPACKH(CONJ(Zeta),Zeta);
      Complex *ZetaL0=ZetaL-k;
      for(; k < stop; ++k) {
        Vec Zetak=ZMULT(X,Y,LOAD(ZetaL0+k));
#else      
    for(unsigned int k=0; k < m;) {
#endif      
#ifdef __SSE2__      
      STORE(u+k,ZMULT(Zetak,LOAD(f+k)));
      STORE(v+k,ZMULT(Zetak,LOAD(g+k)));
#if !SHORTTABLE
      ++k;
      Zetak=LOAD(Zeta+k);
#endif      
#else      
      Complex *P=u+k;
      Complex *Q=v+k;
      Complex fk=*(f+k);
      Complex gk=*(g+k);
      P->re=re*fk.re-im*fk.im;
      P->im=im*fk.re+re*fk.im;
      Q->re=re*gk.re-im*gk.im;
      Q->im=im*gk.re+re*gk.im;
      ++k;
      re=Zeta[k].re;     
      im=Zeta[k].im;     
#endif      
#if SHORTTABLE      
    }  
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
    const Complex ninv2(ninv,ninv);
    Vec Ninv=LOAD(&ninv2);
    Zetak=LOAD(&one);
#else    
    re=1.0;
    im=0.0;
#endif    
#if SHORTTABLE      
    for(unsigned int a=0, k=0; k < m; ++a) {
      unsigned int stop=min(k+s,m);
      Vec Zeta=LOAD(ZetaH+a);
      Vec X=UNPACKL(Zeta,Zeta);
      Vec Y=UNPACKH(CONJ(Zeta),Zeta);
      Complex *ZetaL0=ZetaL-k;
      for(; k < stop; ++k) {
        Vec Zetak=ZMULT(X,Y,LOAD(ZetaL0+k));
#else
    for(unsigned int k=0; k < m;) {
#endif      
#ifdef __SSE2__
      STORE(f+k,Ninv*(ZMULTC(Zetak,LOAD(u+k))+LOAD(f+k)));
#if !SHORTTABLE      
      ++k;
      Zetak=LOAD(Zeta+k);
#endif      
#else      
      Complex *p=f+k;
      Complex fk=*p;
      Complex fkm=*(u+k);
      p->re=ninv*fk.re+re*fkm.re-im*fkm.im;
      p->im=ninv*fk.im+im*fkm.re+re*fkm.im;
      re=Zeta[k];
      im=Zeta[k];
#endif      
    }
#if SHORTTABLE      
    }  
#endif      
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
    
  void padBackwards(Complex *f) {
    unsigned int n2=n/2;
    for(unsigned int i=m; i <= n2; i++) f[i]=0.0;
    cr->fft(f);
  }
  
// Compute f (*) g, where f and g contain the m non-negative Fourier
// components of real functions. Dealiasing is internally implemented via
// explicit zero-padding to size n >= 3*m.
//
// The (distinct) input arrays f and g must each be allocated to size n/2+1
// (contents not preserved). The output is returned in the first m elements
// of f.
  void convolve(Complex *f, Complex *g) {
    padBackwards(f);
    padBackwards(g);
	
    double *F=(double *) f;
    double *G=(double *) g;
    
    double ninv=1.0/n;
    for(unsigned int k=0; k < n; ++k)
      F[k] *= G[k]*ninv;
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
    f[m1]=f0;
    Complex gc=g[c];
    Complex gmk=g[m1];
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
    double fmkre=fmk.re;
    double fmkim=fmk.im;
    double gmkre=gmk.re;
    double gmkim=gmk.im;
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
    for(unsigned int i=0; i < m; ++i)
      V[i] *= U[i];
    rco->fft(v,U); // v is now free

    // r=0:
    V=(double *) v;
    cro->fft(f,V);
    double *F=(double *) f;
    cro->fft(g,F);
    for(unsigned int i=0; i < m; ++i)
      V[i] *= F[i];
    rco->fft(V,f);
    unsigned int cm1=c-1;
    Complex overlap0=f[cm1];
    double overlap1=f[c].re;

    // r=1:
    Complex *f1=f+cm1;
    f1[0]=A;
    f1[1]=fc;
    Complex *g1=g+cm1;
    g1[0]=C;
    g1[1]=gc;
    V=(double *) v;
    cro->fft(g1,V);
    G=(double *) g1;
    cro->fft(f1,G);
    for(unsigned int i=0; i < m; ++i)
      G[i] *= V[i];
    rco->fft(G,v);

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
    for(unsigned k=1; k < cm1; ++k) {
      Complex *p=f+k;
      Complex *s=fm-k;
#ifdef __SSE2__      
      Vec F0=LOAD(p)*ninv2;
      Vec F1=ZMULTC(Zetak,LOAD(v+k));
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
    Complex f1k=conj(Zetak0)*v[cm1];
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
         unsigned int stride, Complex *f) : n(2*m), m(m), M(M), stride(stride) {
    double arg=M_PI/m;
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
      im=im*Cos-re*Sin;
      re=temp;
#endif     
    }
  }
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
    : n(3*m), m(m), M(M), stride(stride) {
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
  
  void padBackwards(Complex *f) {
    for(unsigned int i=0; i < mx; ++i) {
      unsigned int nyi=ny*i;
      unsigned int stop=nyi+ny;
      for(unsigned int j=nyi+my; j < stop; ++j)
        f[j]=0.0;
    }
    
    for(unsigned int i=mx; i < nx; ++i) {
      unsigned int nyi=ny*i;
      unsigned int stop=nyi+ny;
      for(unsigned int j=nyi; j < stop; ++j)
        f[j]=0.0;
    }
    
    if(prune) {
      xBackwards->fft(f);
      yBackwards->fft(f);
    } else
      Backwards->fft(f);
  }
  
  void convolve(Complex *f, Complex *g) {
    padBackwards(f);
    padBackwards(g);
    
    unsigned int n=nx*ny;
    double ninv=1.0/n;
#ifdef __SSE2__      
    const Complex ninv2(ninv,ninv);
    Vec Ninv=LOAD(&ninv2);
    for(unsigned int k=0; k < n; ++k)
      STORE(f+k,Ninv*ZMULT(LOAD(f+k),LOAD(g+k)));
#else    
    for(unsigned int k=0; k < n; ++k)
      f[k] *= g[k]*ninv;
#endif    
	
    if(prune) {
      yForwards->fft(f);
      xForwards->fft(f);
    } else
      Forwards->fft(f);
  }
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
  
  void convolve(Complex *h, Complex *f, Complex *g) {
    for(unsigned int i=0; i < mx; ++i) {
      for(unsigned int j=0; j < my; ++j) {
        Complex sum=0.0;
        for(unsigned int k=0; k <= i; ++k)
          for(unsigned int p=0; p <= j; ++p)
            sum += f[k*my+p]*g[(i-k)*my+j-p];
        h[i*my+j]=sum;
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
  bool odd;   // Is nx odd?
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
    odd=nx % 2 == 1;
    if(odd) {
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
  
  void padBackwards(Complex *f) {
    unsigned int nyp=ny/2+1;
    unsigned int nx2=nx/2;
    unsigned int end=nx2-mx;
    for(unsigned int i=0; i <= end; ++i) {
      unsigned int nypi=nyp*i;
      unsigned int stop=nypi+nyp;
      for(unsigned int j=nypi; j < stop; ++j)
        f[j]=0.0;
    }
    
    for(unsigned int i=nx2+mx; i < nx; ++i) {
      unsigned int nypi=nyp*i;
      unsigned int stop=nypi+nyp;
      for(unsigned int j=nypi; j < stop; ++j)
        f[j]=0.0;
    }
    for(unsigned int i=0; i < nx; ++i) {
      unsigned int nypi=nyp*i;
      unsigned int stop=nypi+nyp;
      for(unsigned int j=nypi+my; j < stop; ++j)
        f[j]=0.0;
    }
    
    if(prune) {
      xBackwards->fft(f);
      if(odd) fftw::Shift(f,nx,ny,-1);
      yBackwards->fft(f);
    } else
      Backwards->fft(f);
  }
  
  void convolve(Complex *f, Complex *g) {
    padBackwards(f);
    padBackwards(g);
    
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

// Out-of-place direct 2D Hermitian convolution.
class DirectHConvolution2 {
protected:  
  unsigned int mx,my;
  
public:
  DirectHConvolution2(unsigned int mx, unsigned int my) : mx(mx), my(my) {}
  
  void convolve(Complex *h, Complex *f, Complex *g) {
    unsigned int xorigin=mx-1;
    int xstart=-xorigin;
    int xstop=mx;
    int ystop=my;
    int ystart=1-my;
    for(int kx=xstart; kx < xstop; ++kx) {
      for(int ky=0; ky < ystop; ++ky) {
        Complex sum=0.0;
        for(int px=xstart; px < xstop; ++px) {
          for(int py=ystart; py < ystop; ++py) {
            int qx=kx-px;
            if(qx >= xstart && qx < xstop) {
              int qy=ky-py;
              if(qy >= ystart && qy < ystop) {
                sum += ((py >= 0) ? f[(xorigin+px)*my+py] : 
                        conj(f[(xorigin-px)*my-py])) *
                  ((qy >= 0) ? g[(xorigin+qx)*my+qy] : 
                   conj(g[(xorigin-qx)*my-qy]));
              }
            }
          }
          h[(xorigin+kx)*my+ky]=sum;
        }
      }
    }	
  }
};

// In-place explicitly dealiased 3D complex convolution.
class ExplicitConvolution3 {
protected:
  unsigned int nx,ny,nz;
  unsigned int mx,my,mz;
  unsigned int nxy;
  unsigned int nyz;
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
    nxy=nx*ny;
    nyz=ny*nz;
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
  
  void padBackwards(Complex *f) {
    for(unsigned int i=0; i < mx; ++i) {
      unsigned int nyi=ny*i;
      for(unsigned int j=0; j < my; ++j) {
        unsigned int nyzij=nz*(nyi+j);
        unsigned int stop=nyzij+nz;
        for(unsigned int k=nyzij+mz; k < stop; ++k)
          f[k]=0.0;
      }
    }
    
    for(unsigned int i=mx; i < nx; ++i) {
      unsigned int nyzi=nyz*i;
      for(unsigned int j=0; j < ny; ++j) {
        unsigned int nyzij=nyzi+nz*j;
        unsigned int stop=nyzij+nz;
        for(unsigned int k=nyzij; k < stop; ++k)
          f[k]=0.0;
      }
    }
    
    for(unsigned int i=0; i < nx; ++i) {
      unsigned int nyzi=nyz*i;
      for(unsigned int j=mx; j < ny; ++j) {
        unsigned int nyzij=nyzi+nz*j;
        unsigned int stop=nyzij+nz;
        for(unsigned int k=nyzij; k < stop; ++k)
          f[k]=0.0;
      }
    }

    if(prune) {
      for(unsigned int i=0; i < mx; i++)
        yBackwards->fft(f+i*nyz);
      xBackwards->fft(f);
      zBackwards->fft(f);
    } else
      Backwards->fft(f);
  }
  
  void convolve(Complex *f, Complex *g) {
    padBackwards(f);
    padBackwards(g);
    
    unsigned int n=nxy*nz;
    double ninv=1.0/n;
#ifdef __SSE2__      
    const Complex ninv2(ninv,ninv);
    Vec Ninv=LOAD(&ninv2);
    for(unsigned int k=0; k < n; ++k)
      STORE(f+k,Ninv*ZMULT(LOAD(f+k),LOAD(g+k)));
#else    
    for(unsigned int k=0; k < n; ++k)
      f[k] *= g[k]*ninv;
#endif    
	
    if(prune) {
      zForwards->fft(f);
      xForwards->fft(f);
     for(unsigned int i=0; i < mx; i++)
        yForwards->fft(f+i*nyz);
     } else
    Forwards->fft(f);
  }
};

// In-place implicitly dealiased 2D complex convolution.
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
  
  void convolve(Complex *h, Complex *f, Complex *g) {
    for(unsigned int i=0; i < mx; ++i) {
      for(unsigned int j=0; j < my; ++j) {
        for(unsigned int k=0; k < mz; ++k) {
          Complex sum=0.0;
          for(unsigned int r=0; r <= i; ++r)
            for(unsigned int p=0; p <= j; ++p)
              for(unsigned int q=0; q <= k; ++q)
                sum += f[r*myz+p*mz+q]*g[(i-r)*myz+(j-p)*mz+(k-q)];
          h[i*myz+j*mz+k]=sum;
        }
      }
    }
  }	
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
    
  void padBackwards(Complex *f) {
    unsigned int n2=n/2;
    for(unsigned int i=m; i <= n2; i++) f[i]=0.0;
    cr->fft(f);
  }
  
// Compute the biconvolution of e, f, and g, where e, f, and g contain the
// m non-negative Fourier components of real functions. Dealiasing is
// internally implemented via explicit zero-padding to size n >= 3*m.
//
// The (distinct) input arrays e, f, and g must each be allocated to size n/2+1
// (contents not preserved). The output is returned in the first m elements
// of e.
  void convolve(Complex *e, Complex *f, Complex *g) {
    padBackwards(e);
    padBackwards(f);
    padBackwards(g);
	
    double *E=(double *) e;
    double *F=(double *) f;
    double *G=(double *) g;
    
    double ninv=1.0/n;
    for(unsigned int k=0; k < n; ++k)
      E[k] *= F[k]*G[k]*ninv;

    rc->fft(e);
  }
};

// In-place implicitly dealiased Hermitian biconvolution.
class ImplicitHBiConvolution {
protected:
  unsigned int n,m;
  rcfft1d *rc, *rco;
  crfft1d *cr, *cro;
  double Cos,Sin;
public:  
  
  // u and v are distinct temporary arrays each of size m+1.
  ImplicitHBiConvolution(unsigned int m, Complex *u, Complex *v) : 
    n(4*m), m(m) {
    double arg=2.0*M_PI/n;
    Cos=cos(arg);
    Sin=sin(arg);

    unsigned int twom=2*m;
    
    rc=new rcfft1d(twom,u);
    cr=new crfft1d(twom,u);
    
    double *U=(double *) u;
    rco=new rcfft1d(twom,U,v);
    cro=new crfft1d(twom,v,U);
  }
  
  ~ImplicitHBiConvolution() {
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
                Complex *u, Complex *v, Complex *w) {
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
      STORE(w+k,ZMULT(Zetak,LOAD(h+k)));
      Zetak=ZMULT(CC,SS,Zetak);
#else
      Complex *P=u+k;
      Complex *Q=v+k;
      Complex *R=w+k;
      Complex fk=*(f+k);
      Complex gk=*(g+k);
      Complex hk=*(h+k);
      P->re=re*fk.re-im*fk.im;
      P->im=im*fk.re+re*fk.im;
      Q->re=re*gk.re-im*gk.im;
      Q->im=im*gk.re+re*gk.im;
      R->re=re*hk.re-im*hk.im;
      R->im=im*hk.re+re*hk.im;
      double temp=re*Cos-im*Sin; 
      im=re*Sin+im*Cos;
      re=temp;
#endif      
    }  
    
    // five of eight FFTs are out-of-place
    
    u[m]=0.0;
    cr->fft(u);
    v[m]=0.0;
    cr->fft(v);
    w[m]=0.0;
    cr->fft(w);
    
    double *U=(double *) u;
    double *V=(double *) v;
    double *W=(double *) w;
    
    unsigned int twom=2*m;
    for(unsigned int i=0; i < twom; ++i)
      V[i] *= U[i]*W[i];
    rco->fft(v,U); // v and w are now free

    V=(double *) v;
    f[m]=0.0;
    cro->fft(f,V);
    double *F=(double *) f;
    g[m]=0.0;
    cro->fft(g,F);
    double *G=(double *) g;
    h[m]=0.0;
    cro->fft(h,G);
    for(unsigned int i=0; i < twom; ++i)
      V[i] *= F[i]*G[i];
    rco->fft(V,f);
    
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
      im=im*Cos-re*Sin;
      re=temp;
#endif      
    }
  }
};

// Out-of-place direct 1D Hermitian biconvolution.
class DirectHBiConvolution {
protected:  
  unsigned int m;
public:
  DirectHBiConvolution(unsigned int m) : m(m) {}
  
  void convolve(Complex *h, Complex *e, Complex *f, Complex *g) {
    int stop=m;
    int start=1-m;
    for(int k=0; k < stop; ++k) {
      Complex sum=0.0;
      for(int p=start; p < stop; ++p) {
        Complex E=(p >= 0) ? e[p] : conj(e[-p]);
        for(int q=start; q < stop; ++q) {
          int r=k-p-q;
          if(r >= start && r < stop)
            sum += E*
              ((q >= 0) ? f[q] : conj(f[-q]))*
              ((r >= 0) ? g[r] : conj(g[-r]));
        }
      }
      h[k]=sum;
    }
  }
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
  unsigned int n;
  unsigned int m;
  unsigned int M;
  unsigned int stride;
  unsigned int dist;
  mfft1d *Backwards;
  mfft1d *Forwards;
  double Cos,Sin;
public:  
  fft0bipad(unsigned int m, unsigned int M, unsigned int stride,
            Complex *f) : n(4*m), m(m), M(M), stride(stride) {
    double arg=2.0*M_PI/n;
    Cos=cos(arg);
    Sin=sin(arg);
    
    unsigned int twom=2*m;
    Backwards=new mfft1d(twom,1,M,stride,1,f);
    Forwards=new mfft1d(twom,-1,M,stride,1,f);
  }
  
  ~fft0bipad() {
    delete Forwards;
    delete Backwards;
  }
  
  void backwards(Complex *f, Complex *u) {
#ifdef __SSE2__
    const Complex cc(Cos,Cos);
    const Complex ss(-Sin,Sin);
    Vec CC=LOAD(&cc);
    Vec SS=LOAD(&ss);
    Complex mizeta(Sin,-Cos);
    Vec Zetak=LOAD(&mizeta);
#else
    double re=Sin;
    double im=-Cos;
#endif    
    for(unsigned int i=0; i < M; ++i)
      f[i]=0.0;
    for(unsigned int i=0; i < M; ++i)
      u[i]=0.0;
    
    unsigned int stop=2*m*stride;
    for(unsigned int k=stride; k < stop; k += stride) {
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
    const Complex ninv1(Sin*ninv,Cos*ninv);
    const Complex ninv2(ninv,ninv);
    Vec CC=LOAD(&cc);
    Vec SS=-LOAD(&ss);
    Vec Zetak=LOAD(&ninv1);
    Vec Ninv2=LOAD(&ninv2);
#else
    double re=Sin*ninv;
    double im=Cos*ninv;
#endif    
    unsigned int stop=2*m*stride;
    for(unsigned int k=stride; k < stop; k += stride) {
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
      im=im*Cos-re*Sin;
      re=temp;
#endif     
    }
  }
};
  
// In-place explicitly dealiased 2D Hermitian biconvolution.
class ExplicitHBiConvolution2 {
protected:
  unsigned int nx,ny;
  unsigned int mx,my;
  bool prune; // Skip Fourier transforming rows containing all zeroes?
  bool odd;   // Is nx odd?
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
    odd=nx % 2 == 1;
    if(odd) {
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
  
  void padBackwards(Complex *f) {
    unsigned int nyp=ny/2+1;
    unsigned int nx2=nx/2;
    unsigned int end=nx2-mx;
    for(unsigned int i=0; i <= end; ++i) {
      unsigned int nypi=nyp*i;
      unsigned int stop=nypi+nyp;
      for(unsigned int j=nypi; j < stop; ++j)
        f[j]=0.0;
    }
    
    for(unsigned int i=nx2+mx; i < nx; ++i) {
      unsigned int nypi=nyp*i;
      unsigned int stop=nypi+nyp;
      for(unsigned int j=nypi; j < stop; ++j)
        f[j]=0.0;
    }
    for(unsigned int i=0; i < nx; ++i) {
      unsigned int nypi=nyp*i;
      unsigned int stop=nypi+nyp;
      for(unsigned int j=nypi+my; j < stop; ++j)
        f[j]=0.0;
    }

    if(prune) {
      xBackwards->fft(f);
      if(odd) fftw::Shift(f,nx,ny,-1);
      yBackwards->fft(f);
    } else
      return Backwards->fft(f);
  }
  
  void convolve(Complex *e, Complex *f, Complex *g) {
    padBackwards(e);
    padBackwards(f);
    padBackwards(g);
    
    double *E=(double *) e;
    double *F=(double *) f;
    double *G=(double *) g;
    
    double ninv=1.0/(nx*ny);
    unsigned int nyp=ny/2+1;
    unsigned int nyp2=2*nyp;

    for(unsigned int i=0; i < nx; ++i) {
      unsigned int nyp2i=nyp2*i;
      unsigned int stop=nyp2i+ny;
      for(unsigned int j=nyp2i; j < stop; ++j)
        E[j] *= F[j]*G[j]*ninv;
    }
	
    if(prune) {
      yForwards->fft(e);
      if(odd) fftw::Shift(e,nx,ny,1);
      xForwards->fft(e);
    } else
      Forwards->fft(e);
  }
};

// In-place implicitly dealiased 2D Hermitian biconvolution.
class ImplicitHBiConvolution2 {
protected:
  unsigned int nx,ny;
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
    xfftpad->backwards(f,u2);
    xfftpad->backwards(g,v2);
    xfftpad->backwards(h,w2);

    unsigned int my1=my+1;
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
  
  void convolve(Complex *h, Complex *e, Complex *f, Complex *g) {
    unsigned int xorigin=mx-1;
    int xstart=-xorigin;
    int xstop=mx;
    int ystop=my;
    int ystart=1-my;
    for(int kx=xstart; kx < xstop; ++kx) {
      for(int ky=0; ky < ystop; ++ky) {
        Complex sum=0.0;
        for(int px=xstart; px < xstop; ++px) {
          for(int py=ystart; py < ystop; ++py) {
            Complex E=(py >= 0) ? e[(xorigin+px)*my+py] : 
              conj(e[(xorigin-px)*my-py]);
            for(int qx=xstart; qx < xstop; ++qx) {
              for(int qy=ystart; qy < ystop; ++qy) {
                int rx=kx-px-qx;
                if(rx >= xstart && rx < xstop) {
                  int ry=ky-py-qy;
                  if(ry >= ystart && ry < ystop) {
                    sum += E *
                      ((qy >= 0) ? f[(xorigin+qx)*my+qy] : 
                       conj(f[(xorigin-qx)*my-qy])) *
                      ((ry >= 0) ? g[(xorigin+rx)*my+ry] : 
                       conj(g[(xorigin-rx)*my-ry]));
                  }
                }
              }
            }
          }
          h[(xorigin+kx)*my+ky]=sum;
        }
      }
    }	
  }
};

#endif
