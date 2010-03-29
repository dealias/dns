#include "Complex.h"
#include "convolution.h"

#ifdef __SSE2__
const union uvec sse2_pm = {
  { 0x00000000,0x00000000,0x00000000,0x80000000 }
};

const union uvec sse2_mm = {
  { 0x00000000,0x80000000,0x00000000,0x80000000 }
};
#endif

const double sqrt3=sqrt(3.0);
const double hsqrt3=0.5*sqrt3;

const Complex hSqrt3(hsqrt3,hsqrt3);
const Complex mhsqrt3(-hsqrt3,-hsqrt3);
const Complex mhalf(-0.5,-0.5);
const Complex zeta3(-0.5,hsqrt3);

unsigned int BuildZeta(unsigned int n, unsigned int m,
                       Complex *&ZetaH, Complex *&ZetaL)
{
  unsigned int s=(int) sqrt(m);
  unsigned int t=m/s;
  if(s*t < m) ++t;
  static const double twopi=2.0*M_PI;
  double arg=twopi/n;
  ZetaH=ComplexAlign(t);
  for(unsigned int a=0; a < t; ++a) {
    double theta=s*a*arg;
    ZetaH[a]=Complex(cos(theta),sin(theta));
  }
  ZetaL=ComplexAlign(s);
  for(unsigned int b=0; b < s; ++b) {
    double theta=b*arg;
    ZetaL[b]=Complex(cos(theta),sin(theta));
  }
  return s;
}

void ExplicitConvolution::padBackwards(Complex *f)
{
  for(unsigned int k=m; k < n; ++k) f[k]=0.0;
  Backwards->fft(f);
}
  
void ExplicitConvolution::convolve(Complex *f, Complex *g)
{
  padBackwards(f);
  padBackwards(g);
      
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

void ImplicitConvolution::preconvolve0(Complex *f, Complex *g, Complex *u,
                                       Complex *v)
{
  // all six FFTs are out-of-place
    
  Backwards->fft(f,u);
  Backwards->fft(g,v);
  
  for(unsigned int k=0; k < m; ++k) {
    Complex *p=u+k;
#ifdef __SSE2__      
    STORE(p,ZMULT(LOAD(p),LOAD(v+k)));
#else      
    Complex uk=*p;
    Complex vk=*(v+k);
    p->re=uk.re*vk.re-uk.im*vk.im;
    p->im=uk.re*vk.im+uk.im*vk.re;
#endif      
  }

  for(unsigned int a=0, k=0; k < m; ++a) {
    unsigned int stop=min(k+s,m);
    Complex *ZetaL0=ZetaL-k;
#ifdef __SSE2__      
    Vec Zeta=LOAD(ZetaH+a);
    Vec X=UNPACKL(Zeta,Zeta);
    Vec Y=UNPACKH(CONJ(Zeta),Zeta);
    for(; k < stop; ++k) {
      Vec Zetak=ZMULT(X,Y,LOAD(ZetaL0+k));
      Complex *fk=f+k;
      Complex *gk=g+k;
      STORE(fk,ZMULT(Zetak,LOAD(fk)));
      STORE(gk,ZMULT(Zetak,LOAD(gk)));
    }
#else
    Complex *p=ZetaH+a;
    double Hre=p->re;
    double Him=p->im;
    for(; k < stop; ++k) {
      Complex *P=f+k;
      Complex *Q=g+k;
      Complex fk=*P;
      Complex gk=*Q;
      Complex L=*(ZetaL0+k);
      double Re=Hre*L.re-Him*L.im;
      double Im=Hre*L.im+Him*L.re;
      P->re=Re*fk.re-Im*fk.im;
      P->im=Im*fk.re+Re*fk.im;
      Q->re=Re*gk.re-Im*gk.im;
      Q->im=Im*gk.re+Re*gk.im;
    }
#endif      
  }
    
  Backwards->fft(f,v);
  Backwards->fft(g,f);
}

void ImplicitConvolution::preconvolve(Complex *f, Complex *g, Complex *u,
                                      Complex *v)
{
  preconvolve0(f,g,u,v);
  for(unsigned int k=0; k < m; ++k) {
#ifdef __SSE2__      
    Complex *vk=v+k;
    STORE(vk,ZMULT(LOAD(vk),LOAD(f+k)));
#else      
    Complex *p=v+k;
    Complex vk=*p;
    Complex fk=*(f+k);
    p->re=vk.re*fk.re-vk.im*fk.im;
    p->im=vk.re*fk.im+vk.im*fk.re;
#endif      
  }
}

void ImplicitConvolution::preconvolvefg(Complex *f, Complex *g, Complex *u,
                                        Complex *v)
{
  preconvolve0(f,g,u,v);
  for(unsigned int k=0; k < m; ++k) {
#ifdef __SSE2__      
    Complex *fk=f+k;
    STORE(fk,ZMULT(LOAD(v+k),LOAD(fk)));
#else      
    Complex *p=f+k;
    Complex fk=*p;
    Complex vk=*(v+k);
    p->re=vk.re*fk.re-vk.im*fk.im;
    p->im=vk.re*fk.im+vk.im*fk.re;
#endif      
  }
  for(unsigned int k=0; k < m; ++k)
    g[k]=u[k];
}

void ImplicitConvolution::postconvolve0(Complex *f, Complex *g, Complex *u,
                                        Complex *v) 
{
  double ninv=0.5/m;
#ifdef __SSE2__      
  const Complex ninv2(ninv,ninv);
  Vec Ninv=LOAD(&ninv2);
#endif    
  for(unsigned int a=0, k=0; k < m; ++a) {
    unsigned int stop=min(k+s,m);
    Complex *ZetaL0=ZetaL-k;
#ifdef __SSE2__
    Vec Zeta=Ninv*LOAD(ZetaH+a);
    Vec X=UNPACKL(Zeta,Zeta);
    Vec Y=UNPACKH(CONJ(Zeta),Zeta);
    for(; k < stop; ++k) {
      Vec Zetak=ZMULT(X,Y,LOAD(ZetaL0+k));
      Complex *fk=f+k;
      STORE(fk,ZMULTC(Zetak,LOAD(u+k))+Ninv*LOAD(fk));
    }
#else      
    Complex *p=ZetaH+a;
    double Hre=ninv*p->re;
    double Him=ninv*p->im;
    for(; k < stop; ++k) {
      Complex *p=f+k;
      Complex fk=*p;
      Complex fkm=*(u+k);
      Complex L=*(ZetaL0+k);
      double Re=Hre*L.re-Him*L.im;
      double Im=Him*L.re+Hre*L.im;
      p->re=ninv*fk.re+Re*fkm.re+Im*fkm.im;
      p->im=ninv*fk.im-Im*fkm.re+Re*fkm.im;
    }
#endif
  }  
}

void ImplicitConvolution::postconvolve(Complex *f, Complex *g, Complex *u,
                                       Complex *v) 
{
  Forwards->fft(u,f);
  Forwards->fft(v,u);
  
  postconvolve0(f,g,u,v);
}

void ImplicitConvolution::postconvolvefg(Complex *f, Complex *g, Complex *u,
                                         Complex *v) 
{
  Forwards->fft(f,u);
  Forwards->fft(g,f);
  
  postconvolve0(f,g,u,v);
}

void ImplicitConvolution::convolve(Complex *f, Complex *g, Complex *u,
                                   Complex *v) 
{
  preconvolve(f,g,u,v);
  postconvolve(f,g,u,v);
}

void DirectConvolution::convolve(Complex *h, Complex *f, Complex *g)
{
  for(unsigned int i=0; i < m; ++i) {
    Complex sum=0.0;
    for(unsigned int j=0; j <= i; ++j) sum += f[j]*g[i-j];
    h[i]=sum;
  }
}	

void ExplicitHConvolution::padBackwards(Complex *f)
{
  unsigned int n2=n/2;
  for(unsigned int i=m; i <= n2; ++i) f[i]=0.0;
  cr->fft(f);
}
  
void ExplicitHConvolution::convolve(Complex *f, Complex *g)
{
  padBackwards(f);
  padBackwards(g);
	
  double *F=(double *) f;
  double *G=(double *) g;
    
  double ninv=1.0/n;
  for(unsigned int k=0; k < n; ++k)
    F[k] *= G[k]*ninv;
  rc->fft(f);
}

void ImplicitHConvolution::convolve(Complex *f, Complex *g, Complex *u,
                                    Complex *v)
{
  double f0=f[0].re;
  double g0=g[0].re;

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
  unsigned int stop=s;
  Complex *ZetaL0=ZetaL;
    
#ifdef __SSE2__      
  Vec Fmk=LOAD(&fmk);
  Vec Gmk=LOAD(&gmk);
  Vec Mhalf=LOAD(&mhalf);
  Vec HSqrt3=LOAD(&hSqrt3);
  for(unsigned int a=0, k=1; k < c; ++a) {
    Vec Zeta=LOAD(ZetaH+a);
    Vec X=UNPACKL(Zeta,Zeta);
    Vec Y=UNPACKH(CONJ(Zeta),Zeta);
    for(; k < stop; ++k) {
      Vec Zetak=ZMULT(X,Y,LOAD(ZetaL0+k));
      Complex *p=f+k;
      Complex *q=g+k;
      Vec A=LOAD(p);
      Vec B=LOAD(q);
      Vec C=Fmk*Mhalf+CONJ(A);
      Vec D=Gmk*Mhalf+CONJ(B);
      STORE(p,A+CONJ(Fmk));
      STORE(q,B+CONJ(Gmk));
      Fmk *= HSqrt3;
      Gmk *= HSqrt3;
      A=ZMULTC(Zetak,UNPACKL(C,Fmk));
      B=ZMULTIC(Zetak,UNPACKH(C,Fmk));
      C=ZMULTC(Zetak,UNPACKL(D,Gmk));
      D=ZMULTIC(Zetak,UNPACKH(D,Gmk));
      STORE(u+k,A-B);
      STORE(v+k,C-D);
      p=f+m1-k;
      Fmk=LOAD(p);
      STORE(p,A+B);
      q=g+m1-k;
      Gmk=LOAD(q);
      STORE(q,C+D);
    }
    stop=min(k+s,c);
    ZetaL0=ZetaL-k;
  }
#else
  for(unsigned int a=0, k=1; k < c; ++a) {
    double fmkre=fmk.re;
    double fmkim=fmk.im;
    double gmkre=gmk.re;
    double gmkim=gmk.im;
    Complex *p=ZetaH+a;
    double Hre=p->re;
    double Him=-p->im;
    for(; k < stop; ++k) {
      Complex *p=f+k;
      Complex *q=g+k;
      Complex L=*(ZetaL0+k);
      double Re=Hre*L.re+Him*L.im;
      double Im=Him*L.re-Hre*L.im;
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
    }
    stop=min(k+s,c);
    ZetaL0=ZetaL-k;
  }
#endif      
  
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
  double *G=(double *) g1;
  cro->fft(f1,G);
  for(unsigned int i=0; i < m; ++i)
    G[i] *= V[i];
  rco->fft(G,v);

  double ninv=1.0/(3.0*m);
  f[0]=(f[0].re+v[0].re+u[0].re)*ninv;
  Complex *fm=f+m;
  stop=s;
  ZetaL0=ZetaL;
#ifdef __SSE2__      
  const Complex ninv2(ninv,ninv);
  Vec Ninv=LOAD(&ninv2);
  Fmk=LOAD(&fmk);
  Gmk=LOAD(&gmk);
  for(unsigned int a=0, k=1; k < cm1; ++a) {
    Vec Zeta=Ninv*LOAD(ZetaH+a);
    Vec X=UNPACKL(Zeta,Zeta);
    Vec Y=UNPACKH(CONJ(Zeta),Zeta);
    for(; k < stop; ++k) {
      Complex *p=f+k;
      Complex *s=fm-k;
      Vec Zetak=ZMULT(X,Y,LOAD(ZetaL0+k));
      Vec F0=LOAD(p)*Ninv;
      Vec F1=ZMULTC(Zetak,LOAD(v+k));
      Vec F2=ZMULT(Zetak,LOAD(u+k));
      Vec S=F1+F2;
      STORE(p,F0+S);
      STORE(s,CONJ(F0+Mhalf*S)-HSqrt3*FLIP(F1-F2));
    }
    stop=min(k+s,cm1);
    ZetaL0=ZetaL-k;
  }
#else
  for(unsigned int a=0, k=1; k < cm1; ++a) {
    Complex *p=ZetaH+a;
    double Hre=ninv*p->re;
    double Him=ninv*p->im;
    for(; k < stop; ++k) {
      Complex *p=f+k;
      Complex *s=fm-k;
      Complex q=v[k];
      Complex r=u[k];
      Complex L=*(ZetaL0+k);
      double Re=Hre*L.re-Him*L.im;
      double Im=Him*L.re+Hre*L.im;
      double f0re=p->re*ninv;
      double f0im=p->im*ninv;
      double f1re=Re*q.re+Im*q.im;
      double f2re=Re*r.re-Im*r.im;
      double sre=f1re+f2re;
      double f1im=Re*q.im-Im*q.re;
      double f2im=Re*r.im+Im*r.re;
      double sim=f1im+f2im;
      p->re=f0re+sre;
      p->im=f0im+sim;
      s->re=f0re-0.5*sre-hsqrt3*(f1im-f2im);
      s->im=-f0im+0.5*sim-hsqrt3*(f1re-f2re);
    }
    stop=min(k+s,cm1);
    ZetaL0=ZetaL-k;
  }
#endif    
    
  unsigned int a=(c-1)/s;
  Complex Zetak0=ninv*ZetaH[a]*ZetaL[c-1-s*a];
  Complex f0k=overlap0*ninv;
  Complex f1k=conj(Zetak0)*v[cm1];
  Complex f2k=Zetak0*u[cm1];
  f[cm1]=f0k+f1k+f2k;
  if(c > 1) f[c+1]=conj(f0k+zeta3*f1k)+zeta3*conj(f2k);
  f[c]=(overlap1-v[c].re*zeta3-u[c].re*conj(zeta3))*ninv;
}

void DirectHConvolution::convolve(Complex *h, Complex *f, Complex *g)
{
  for(unsigned int i=0; i < m; ++i) {
    Complex sum=0.0;
    for(unsigned int j=0; j <= i; ++j) sum += f[j]*g[i-j];
    for(unsigned int j=i+1; j < m; ++j) sum += f[j]*conj(g[j-i]);
    for(unsigned int j=1; j < m-i; ++j) sum += conj(f[j])*g[i+j];
    h[i]=sum;
  }
}	

void fftpad::backwards(Complex *f, Complex *u)
{
  for(unsigned int a=0, k=0; k < m; ++a) {
    unsigned int stop=min(k+s,m);
    Complex *ZetaL0=ZetaL-k;
#ifdef __SSE2__      
    Vec H=LOAD(ZetaH+a);
    for(; k < stop; ++k) {
      Vec Zetak=ZMULT(H,LOAD(ZetaL0+k));
      Vec X=UNPACKL(Zetak,Zetak);
      Vec Y=UNPACKH(CONJ(Zetak),Zetak);
      unsigned int kstride=k*stride;
      Complex *fk=f+kstride;
      Complex *uk=u+kstride;
      for(unsigned int i=0; i < M; ++i)
        STORE(uk+i,ZMULT(X,Y,LOAD(fk+i)));
    }
#else
    Complex H=ZetaH[a];
    for(; k < stop; ++k) {
      Complex L=ZetaL0[k];
      double Re=H.re*L.re-H.im*L.im;
      double Im=H.re*L.im+H.im*L.re;
      unsigned int kstride=k*stride;
      Complex *fk=f+kstride;
      Complex *uk=u+kstride;
      for(unsigned int i=0; i < M; ++i) {
        Complex *p=uk+i;
        Complex fki=*(fk+i);
        p->re=Re*fki.re-Im*fki.im;
        p->im=Im*fki.re+Re*fki.im;
      }
    }
#endif      
  }
    
  Backwards->fft(f);
  Backwards->fft(u);
}
  
void fftpad::forwards(Complex *f, Complex *u)
{
  Forwards->fft(f);
  Forwards->fft(u);

  double ninv=0.5/m;
#ifdef __SSE2__
  const Complex ninv2(ninv,ninv);
  Vec Ninv=LOAD(&ninv2);
#endif    
  for(unsigned int a=0, k=0; k < m; ++a) {
    unsigned int stop=min(k+s,m);
    Complex *ZetaL0=ZetaL-k;
#ifdef __SSE2__      
    Vec H=Ninv*LOAD(ZetaH+a);
    for(; k < stop; ++k) {
      Vec Zetak=ZMULT(H,LOAD(ZetaL0+k));
      Vec X=UNPACKL(Zetak,Zetak);
      Vec Y=UNPACKH(Zetak,CONJ(Zetak));
      unsigned int kstride=k*stride;
      Complex *uk=u+kstride;
      Complex *fk=f+kstride;
      for(unsigned int i=0; i < M; ++i) {
        Complex *p=fk+i;
        STORE(p,LOAD(p)*Ninv+ZMULT(X,Y,LOAD(uk+i)));
      }
    }
#else
    Complex H=ninv*ZetaH[a];
    for(; k < stop; ++k) {
      Complex L=ZetaL0[k];
      double Re=H.re*L.re-H.im*L.im;
      double Im=H.re*L.im+H.im*L.re;
      unsigned int kstride=k*stride;
      Complex *uk=u+kstride;
      Complex *fk=f+kstride;
      for(unsigned int i=0; i < M; ++i) {
        Complex *p=fk+i;
        Complex fki=*p;
        Complex fkm=*(uk+i);
        p->re=ninv*fki.re+Re*fkm.re+Im*fkm.im;
        p->im=ninv*fki.im-Im*fkm.re+Re*fkm.im;
      }
    }
#endif     
  }
}

void fft0pad::backwards(Complex *f, Complex *u)
{
  unsigned int m1=m-1;
  unsigned int m1stride=m1*stride;
  Complex *fm1stride=f+m1stride;
  for(unsigned int i=0; i < M; ++i)
    u[i]=fm1stride[i];
    
#ifdef __SSE2__      
  Vec Mhalf=LOAD(&mhalf);
  Vec Mhsqrt3=LOAD(&mhsqrt3);
#endif    
  unsigned int stop=s;
  Complex *ZetaL0=ZetaL;
  for(unsigned int a=0, k=1; k < m; ++a) {
#ifdef __SSE2__
    Vec H=LOAD(ZetaH+a);
    for(; k < stop; ++k) {
      Vec Zetak=ZMULT(H,LOAD(ZetaL0+k));
      Vec X=UNPACKL(Zetak,Zetak);
      Vec Y=UNPACKH(CONJ(Zetak),Zetak);
      unsigned int kstride=k*stride;
      Complex *uk=u+kstride;
      Complex *fk=f+kstride;
      Complex *fmk=fm1stride+kstride;
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
    }
#else        
    Complex H=ZetaH[a];
    for(; k < stop; ++k) {
      Complex L=ZetaL0[k];
      double Re=H.re*L.re-H.im*L.im;
      double Im=H.re*L.im+H.im*L.re;
      unsigned int kstride=k*stride;
      Complex *uk=u+kstride;
      Complex *fk=f+kstride;
      Complex *fmk=fm1stride+kstride;
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
    }
#endif
    stop=min(k+s,m);
    ZetaL0=ZetaL-k;
  }
    
  Backwards->fft(f);
  Complex *umstride=u+m*stride;
  for(unsigned int i=0; i < M; ++i) {
    umstride[i]=fm1stride[i]; // Store extra value here.
    fm1stride[i]=u[i];
  }
    
  Backwards->fft(fm1stride);
  Backwards->fft(u);
}
  

void fft0pad::forwards(Complex *f, Complex *u)
{
  unsigned int m1stride=(m-1)*stride;
  Complex *fm1stride=f+m1stride;
  Forwards->fft(fm1stride);
  Complex *umstride=u+m*stride;
  for(unsigned int i=0; i < M; ++i) {
    Complex temp=umstride[i];
    umstride[i]=fm1stride[i];
    fm1stride[i]=temp;
  }
    
  Forwards->fft(f);
  Forwards->fft(u);

  double ninv=1.0/(3.0*m);
  for(unsigned int i=0; i < M; ++i)
    umstride[i]=(umstride[i]+f[i]+u[i])*ninv;
#ifdef __SSE2__      
  const Complex ninv2(ninv,ninv);
  Vec Ninv=LOAD(&ninv2);
  Vec Mhalf=LOAD(&mhalf);
  Vec HSqrt3=LOAD(&hSqrt3);
#endif    
  unsigned int stop=s;
  Complex *ZetaL0=ZetaL;
  for(unsigned int a=0, k=1; k < m; ++a) {
#ifdef __SSE2__      
    Vec H=LOAD(ZetaH+a)*Ninv;
    for(; k < stop; ++k) {
      Vec Zetak=ZMULT(H,LOAD(ZetaL0+k));
      Vec X=UNPACKL(Zetak,Zetak);
      Vec Y=UNPACKH(CONJ(Zetak),Zetak);
      unsigned int kstride=k*stride;
      Complex *fk=f+kstride;
      Complex *fm1k=fm1stride+kstride;
      Complex *uk=u+kstride;
      for(unsigned int i=0; i < M; ++i) {
        Complex *p=fk+i;
        Complex *q=fm1k+i;
        Vec F0=LOAD(p)*Ninv;
        Vec F1=ZMULT(X,-Y,LOAD(q));
        Vec F2=ZMULT(X,Y,LOAD(uk+i));
        Vec S=F1+F2;
        STORE(p-stride,F0+Mhalf*S+HSqrt3*ZMULTI(F1-F2));
        STORE(q,F0+S);
      }
    }
#else
    Complex H=ninv*ZetaH[a];
    for(; k < stop; ++k) {
      Complex L=ZetaL0[k];
      double Re=H.re*L.re-H.im*L.im;
      double Im=H.re*L.im+H.im*L.re;
      unsigned int kstride=k*stride;
      Complex *fk=f+kstride;
      Complex *fm1k=fm1stride+kstride;
      Complex *uk=u+kstride;
      for(unsigned int i=0; i < M; ++i) {
        Complex *p=fk+i;
        Complex *q=fm1k+i;
        Complex z=*q;
        Complex r=uk[i];
        double f0re=p->re*ninv;
        double f0im=p->im*ninv;
        double f1re=Re*z.re+Im*z.im;
        double f1im=Re*z.im-Im*z.re;
        double f2re=Re*r.re-Im*r.im;
        double f2im=Re*r.im+Im*r.re;
        double sre=f1re+f2re;
        double sim=f1im+f2im;
        p -= stride;
        p->re=f0re-0.5*sre-hsqrt3*(f1im-f2im);
        p->im=f0im-0.5*sim+hsqrt3*(f1re-f2re);
        q->re=f0re+sre;
        q->im=f0im+sim;
      }
    }
#endif      
    stop=min(k+s,m);
    ZetaL0=ZetaL-k;
  }
  for(unsigned int i=0; i < M; ++i)
    fm1stride[i]=umstride[i];
}

void ExplicitConvolution2::padBackwards(Complex *f)
{
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
  
void ExplicitConvolution2::convolve(Complex *f, Complex *g)
{
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

void DirectConvolution2::convolve(Complex *h, Complex *f, Complex *g)
{
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

void ExplicitHConvolution2::padBackwards(Complex *f)
{
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
    if(nx % 2) fftw::Shift(f,nx,ny,-1);
    yBackwards->fft(f);
  } else
    Backwards->fft(f);
}

void ExplicitHConvolution2::convolve(Complex *f, Complex *g)
{
  unsigned int xorigin=nx/2;
  unsigned int nyp=ny/2+1;
    
  HermitianSymmetrizeX(mx,nyp,xorigin,f);
  HermitianSymmetrizeX(mx,nyp,xorigin,g);
    
  padBackwards(f);
  padBackwards(g);
    
  double *F=(double *) f;
  double *G=(double *) g;
    
  double ninv=1.0/(nx*ny);
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

void DirectHConvolution2::convolve(Complex *h, Complex *f, Complex *g)
{
  unsigned int xorigin=mx-1;
    
  HermitianSymmetrizeX(mx,my,xorigin,f);
  HermitianSymmetrizeX(mx,my,xorigin,g);
    
  int xstart=-xorigin;
  int ystart=1-my;
  int xstop=mx;
  int ystop=my;
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

void ExplicitConvolution3::padBackwards(Complex *f)
{
  for(unsigned int i=0; i < mx; ++i) {
    unsigned int nyi=ny*i;
    for(unsigned int j=0; j < my; ++j) {
      unsigned int nyzij=nz*(nyi+j);
      unsigned int stop=nyzij+nz;
      for(unsigned int k=nyzij+mz; k < stop; ++k)
        f[k]=0.0;
    }
  }
    
  unsigned int nyz=ny*nz;
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
    for(unsigned int i=0; i < mx; ++i)
      yBackwards->fft(f+i*nyz);
    xBackwards->fft(f);
    zBackwards->fft(f);
  } else
    Backwards->fft(f);
}
  
void ExplicitConvolution3::convolve(Complex *f, Complex *g)
{
  padBackwards(f);
  padBackwards(g);
    
  unsigned int n=nx*ny*nz;
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
    unsigned int nyz=ny*nz;
    for(unsigned int i=0; i < mx; ++i)
      yForwards->fft(f+i*nyz);
  } else
    Forwards->fft(f);
}

void DirectConvolution3::convolve(Complex *h, Complex *f, Complex *g)
{
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

void DirectHConvolution3::convolve(Complex *h, Complex *f, Complex *g)
{
  unsigned int xorigin=mx-1;
  unsigned int yorigin=my-1;
  unsigned int ny=2*my-1;
    
  HermitianSymmetrizeXY(mx,my,mz,ny,xorigin,yorigin,f);
  HermitianSymmetrizeXY(mx,my,mz,ny,xorigin,yorigin,g);
    
  int xstart=-xorigin;
  int ystart=-yorigin;
  int zstart=1-mz;
  int xstop=mx;
  int ystop=my;
  int zstop=mz;
  for(int kx=xstart; kx < xstop; ++kx) {
    for(int ky=ystart; ky < ystop; ++ky) {
      for(int kz=0; kz < zstop; ++kz) {
        Complex sum=0.0;
        for(int px=xstart; px < xstop; ++px) {
          for(int py=ystart; py < ystop; ++py) {
            for(int pz=zstart; pz < zstop; ++pz) {
              int qx=kx-px;
              if(qx >= xstart && qx < xstop) {
                int qy=ky-py;
                if(qy >= ystart && qy < ystop) {
                  int qz=kz-pz;
                  if(qz >= zstart && qz < zstop) {
                    sum += ((pz >= 0) ? 
                            f[((xorigin+px)*ny+yorigin+py)*mz+pz] : 
                            conj(f[((xorigin-px)*ny+yorigin-py)*mz-pz])) *
                      ((qz >= 0) ? g[((xorigin+qx)*ny+yorigin+qy)*mz+qz] :    
                       conj(g[((xorigin-qx)*ny+yorigin-qy)*mz-qz]));
                  }
                }
              }
            }
          }
        }
        h[((xorigin+kx)*ny+yorigin+ky)*mz+kz]=sum;
      }
    }
  }	
}

void ExplicitHBiConvolution::padBackwards(Complex *f)
{
    unsigned int n2=n/2;
    for(unsigned int i=m; i <= n2; ++i) f[i]=0.0;
    cr->fft(f);
}

void ExplicitHBiConvolution::convolve(Complex *f, Complex *g, Complex *h)
{
  padBackwards(f);
  padBackwards(g);
  padBackwards(h);
	
  double *F=(double *) f;
  double *G=(double *) g;
  double *H=(double *) h;
    
  double ninv=1.0/n;
  for(unsigned int k=0; k < n; ++k)
    F[k] *= G[k]*H[k]*ninv;

  rc->fft(f);
}

void ImplicitHBiConvolution::convolve(Complex *f, Complex *g, Complex *h,
                                      Complex *u, Complex *v, Complex *w)
{
  for(unsigned int a=0, k=0; k < m; ++a) {
    unsigned int stop=min(k+s,m);
    Complex *ZetaL0=ZetaL-k;
#ifdef __SSE2__      
    Vec Zeta=LOAD(ZetaH+a);
    Vec X=UNPACKL(Zeta,Zeta);
    Vec Y=UNPACKH(CONJ(Zeta),Zeta);
    for(; k < stop; ++k) {
      Vec Zetak=ZMULT(X,Y,LOAD(ZetaL0+k));
      STORE(u+k,ZMULT(Zetak,LOAD(f+k)));
      STORE(v+k,ZMULT(Zetak,LOAD(g+k)));
      STORE(w+k,ZMULT(Zetak,LOAD(h+k)));
    }
#else
    Complex *p=ZetaH+a;
    double Hre=p->re;
    double Him=p->im;
    for(; k < stop; ++k) {
      Complex *P=u+k;
      Complex *Q=v+k;
      Complex *R=w+k;
      Complex fk=*(f+k);
      Complex gk=*(g+k);
      Complex hk=*(h+k);
      Complex L=*(ZetaL0+k);
      double Re=Hre*L.re-Him*L.im;
      double Im=Hre*L.im+Him*L.re;
      P->re=Re*fk.re-Im*fk.im;
      P->im=Im*fk.re+Re*fk.im;
      Q->re=Re*gk.re-Im*gk.im;
      Q->im=Im*gk.re+Re*gk.im;
      R->re=Re*hk.re-Im*hk.im;
      R->im=Im*hk.re+Re*hk.im;
    }
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
    
  double ninv=0.25/m;
#ifdef __SSE2__      
  const Complex ninv2(ninv,ninv);
  Vec Ninv=LOAD(&ninv2);
#endif    
  for(unsigned int a=0, k=0; k < m; ++a) {
    unsigned int stop=min(k+s,m);
    Complex *ZetaL0=ZetaL-k;
#ifdef __SSE2__
    Vec Zeta=Ninv*LOAD(ZetaH+a);
    Vec X=UNPACKL(Zeta,Zeta);
    Vec Y=UNPACKH(CONJ(Zeta),Zeta);
    for(; k < stop; ++k) {
      Vec Zetak=ZMULT(X,Y,LOAD(ZetaL0+k));
      STORE(f+k,ZMULTC(Zetak,LOAD(u+k))+Ninv*LOAD(f+k));
    }
#else      
    Complex *p=ZetaH+a;
    double Hre=ninv*p->re;
    double Him=ninv*p->im;
    for(; k < stop; ++k) {
      Complex *p=f+k;
      Complex fk=*p;
      Complex fkm=*(u+k);
      Complex L=*(ZetaL0+k);
      double Re=Hre*L.re-Him*L.im;
      double Im=Him*L.re+Hre*L.im;
      p->re=ninv*fk.re+Re*fkm.re+Im*fkm.im;
      p->im=ninv*fk.im-Im*fkm.re+Re*fkm.im;
    }
#endif      
  }
}

void DirectHBiConvolution::convolve(Complex *h, Complex *e, Complex *f,
                                    Complex *g)
{
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

void fft0bipad::backwards(Complex *f, Complex *u)
{
  for(unsigned int i=0; i < M; ++i)
    f[i]=0.0;
  for(unsigned int i=0; i < M; ++i)
    u[i]=0.0;
    
  unsigned int twom=2*m;
  unsigned int stop=s;
  Complex *ZetaL0=ZetaL;
  for(unsigned int a=0, k=1; k < twom; ++a) {
#ifdef __SSE2__      
    Vec H=-LOAD(ZetaH+a);
    for(; k < stop; ++k) {
      Vec Zetak=ZMULT(H,LOAD(ZetaL0+k));
      Vec X=UNPACKL(Zetak,Zetak);
      Vec Y=UNPACKH(CONJ(Zetak),Zetak);
      unsigned int kstride=k*stride;
      Complex *fk=f+kstride;
      Complex *uk=u+kstride;
      for(unsigned int i=0; i < M; ++i)
        STORE(uk+i,ZMULTI(X,Y,LOAD(fk+i)));
    }
#else
    Complex H=ZetaH[a];
    for(; k < stop; ++k) {
      Complex L=ZetaL0[k];
      double Re=H.im*L.re+H.re*L.im;
      double Im=H.im*L.im-H.re*L.re;
      unsigned int kstride=k*stride;
      Complex *fk=f+kstride;
      Complex *uk=u+kstride;
      for(unsigned int i=0; i < M; ++i) {
        Complex *p=uk+i;
        Complex fki=*(fk+i);
        p->re=Re*fki.re-Im*fki.im;
        p->im=Im*fki.re+Re*fki.im;
      }
    }
#endif      
    stop=min(k+s,twom);
    ZetaL0=ZetaL-k;
  }
    
  Backwards->fft(f);
  Backwards->fft(u);
}

void fft0bipad::forwards(Complex *f, Complex *u)
{
  Forwards->fft(f);
  Forwards->fft(u);

  double ninv=0.25/m;
#ifdef __SSE2__
  const Complex ninv2(ninv,ninv);
  Vec Ninv=LOAD(&ninv2);
#endif    
  unsigned int twom=2*m;
  unsigned int stop=s;
  Complex *ZetaL0=ZetaL;
  for(unsigned int a=0, k=1; k < twom; ++a) {
#ifdef __SSE2__      
    Vec H=Ninv*LOAD(ZetaH+a);
    for(; k < stop; ++k) {
      Vec Zetak=ZMULT(H,LOAD(ZetaL0+k));
      Vec X=UNPACKL(Zetak,Zetak);
      Vec Y=UNPACKH(Zetak,CONJ(Zetak));
      unsigned int kstride=k*stride;
      Complex *uk=u+kstride;
      Complex *fk=f+kstride;
      for(unsigned int i=0; i < M; ++i) {
        Complex *p=fk+i;
        STORE(p,LOAD(p)*Ninv+ZMULTI(X,Y,LOAD(uk+i)));
      }
    }
#else        
    Complex H=ninv*ZetaH[a];
    for(; k < stop; ++k) {
      Complex L=ZetaL0[k];
      double Re=H.im*L.re+H.re*L.im;
      double Im=H.re*L.re-H.im*L.im;
      unsigned int kstride=k*stride;
      Complex *uk=u+kstride;
      Complex *fk=f+kstride;
      for(unsigned int i=0; i < M; ++i) {
        Complex *p=fk+i;
        Complex fki=*p;
        Complex fkm=*(uk+i);
        p->re=ninv*fki.re+Re*fkm.re-Im*fkm.im;
        p->im=ninv*fki.im+Im*fkm.re+Re*fkm.im;
      }
    }
#endif     
    stop=min(k+s,twom);
    ZetaL0=ZetaL-k;
  }
}

void ExplicitHBiConvolution2::padBackwards(Complex *f)
{
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
    if(nx % 2) fftw::Shift(f,nx,ny,-1);
    yBackwards->fft(f);
  } else
    return Backwards->fft(f);
}

void ExplicitHBiConvolution2::convolve(Complex *f, Complex *g, Complex *h)
{
  unsigned int xorigin=nx/2;
  unsigned int nyp=ny/2+1;
    
  HermitianSymmetrizeX(mx,nyp,xorigin,f);
  HermitianSymmetrizeX(mx,nyp,xorigin,g);
  HermitianSymmetrizeX(mx,nyp,xorigin,h);
    
  padBackwards(f);
  padBackwards(g);
  padBackwards(h);
    
  double *F=(double *) f;
  double *G=(double *) g;
  double *H=(double *) h;
    
  double ninv=1.0/(nx*ny);
  unsigned int nyp2=2*nyp;

  for(unsigned int i=0; i < nx; ++i) {
    unsigned int nyp2i=nyp2*i;
    unsigned int stop=nyp2i+ny;
    for(unsigned int j=nyp2i; j < stop; ++j)
      F[j] *= G[j]*H[j]*ninv;
  }
	
  if(prune) {
    yForwards->fft(f);
    if(nx % 2) fftw::Shift(f,nx,ny,1);
    xForwards->fft(f);
  } else
    Forwards->fft(f);
}

void DirectHBiConvolution2::convolve(Complex *h, Complex *e, Complex *f,
                                     Complex *g)
{
  HermitianSymmetrizeX(mx,my,mx-1,e);
  HermitianSymmetrizeX(mx,my,mx-1,f);
  HermitianSymmetrizeX(mx,my,mx-1,g);
    
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
