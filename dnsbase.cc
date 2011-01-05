#include "dnsbase.h"

const int DNSBase::xpad=1;
const int DNSBase::ypad=1;

//***** Source routines *****//

static const double sqrt2=sqrt(2.0);

void DNSBase::LinearSource(const vector& wSrc, const vector& w0, double)
{
  w.Set(w0);
  f0.Set(wSrc);

  for(unsigned i=0; i < xorigin; i++) {
    unsigned I=xorigin-i;
    unsigned I2=I*I;
    unsigned im=xorigin+I;
    
    vector f0i=f0[i];
    vector wi=w[i];
    vector f0im=f0[im];
    vector wim=w[im];
    
    // diagonals
    Real nukk=nuk(2*I2);
    f0i[I] -= nukk*wi[I];
    f0im[I] -= nukk*wim[I];
    
    const unsigned stop=I;
    for(unsigned j=1; j < stop; ++j) {
      nukk=nuk(I2+j*j);
      // bottom left
      f0i[j] -= nukk*wi[j];
      // bottom right
      f0im[j] -= nukk*wim[j];
      // top left
      unsigned jm=xorigin-j;
      f0[jm][I] -= nukk*w[jm][I];
      // top right
      unsigned jp=xorigin+j;
      f0[jp][I] -= nukk*w[jp][I];
    }
  }
  
  // xorigin case
  vector f0i=f0[xorigin];
  vector wi=w[xorigin];
  for(unsigned j=1; j < my; ++j)
    f0i[j] -= nuk(j*j)*wi[j];
  
  // bottom right
  for(unsigned i=xorigin; i < Nx; ++i) {
    unsigned I=i-xorigin;
    f0[i][0] -= nuk(I*I)*w[i][0]; 
  }
}

void DNSBase::NonLinearSource(const vector& wSrc, const vector& wY, double)
{
  w.Set(wY);
  f0.Set(wSrc);

  f0(origin)=0.0;
  f1(origin)=0.0;
  g0(origin)=0.0;
  g1(origin)=0.0;
  // TODO: optimize?
  for(unsigned i=0; i < Nx; ++i) {
    Real kx=k0*((int) i-(int) xorigin);
    Real kx2=kx*kx;
    vector wi=w[i];
    vector f0i=f0[i];
    vector f1i=f1[i];
    vector g0i=g0[i];
    vector g1i=g1[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      Real ky=k0*j;
      Complex wij=wi[j];
      Complex kxw=Complex(-kx*wij.im,kx*wij.re);
      Complex kyw=Complex(-ky*wij.im,ky*wij.re);
      f0i[j]=kxw;
      f1i[j]=kyw;
      Real k2inv=1.0/(kx2+ky*ky);
      g0i[j]=k2inv*kyw;
      g1i[j]=-k2inv*kxw;
    }
  }
  F[0]=f0;
  Convolution->convolve(F,G);
  f0(origin)=0.0;
  
#if 0
  Real sum=0.0;
  for(unsigned i=0; i < Nx; ++i) {
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      Complex wij=wi[j];
      sum += (f0[i][j]*conj(wij)).re;
    }
  }
  cout << sum << endl;
#endif

}

void DNSBase::Transfer(const vector2& Src, const vector2& Y)
{
  w.Set(Y[OMEGA]);
  f0.Set(Src[OMEGA]);
  Set(T,Src[TRANSFER]);
  for(unsigned K=0; K < nshells; K++)
    T[K]=Complex(0.0,0.0);

  for(unsigned i=0; i < xorigin; i++) {
    unsigned I=xorigin-i;
    unsigned im=xorigin+xorigin-i;
    unsigned I2=I*I;
    
    vector wi=w[i];
    vector wim=w[im];
    vector fi=f0[i];
    vector fim=f0[im];
    
    // diagonals
    unsigned Sk= spectrum == DISCRETE ?  kval[I][I] : (unsigned)(sqrt2*I-0.5);
    T[Sk].re += realproduct(fi[I],wi[I]) + realproduct(fim[I],wim[I]);
   
    const unsigned stop=I;
    for(unsigned j=1; j < stop; ++j) {
      unsigned jm=xorigin-j;
      unsigned jp=xorigin+j;
      Sk= spectrum == DISCRETE ? kval[I][j] : (unsigned)(sqrt((I2+j*j))-0.5);
      T[Sk].re += 
	realproduct(fi[j],wi[j]) +
	realproduct(fim[j],wim[j]) +
	realproduct(f0[jm][I],w[jm][I]) +
	realproduct(f0[jp][I],w[jp][I]);
    }
  }

  // xorigin case
  vector wi=w[xorigin];
  vector fi=f0[xorigin];
  for(unsigned j=1; j < my; ++j) {
    unsigned Sk= spectrum == DISCRETE ? kval[j][0] : j-1;
    T[Sk].re += realproduct(fi[j],wi[j]);
  }
  // bottom right
  for(unsigned i=xorigin+1; i < Nx; ++i) {
    unsigned Sk= spectrum == DISCRETE ? kval[i-xorigin][0] : i-xorigin-1;
    T[Sk].re += realproduct(f0[i][0],w[i][0]);
  }  
  Forcing->Force(f0,T);
  
  if(casimir)
    CasimirTransfer(Src,Y);
}

void DNSBase::Spectrum(vector& S, const vector& y) 
{
  w.Set(y);
  for(unsigned K=0; K < nshells; K++)    
    S[K]=Complex(0.0,0.0);

  for(unsigned i=0; i < xorigin; i++) {
    unsigned I=xorigin-i;
    unsigned im=xorigin+I;
    unsigned I2=I*I;

    vector wi=w[i];
    vector wim=w[im];

    // diagonals
    unsigned k2=2*I2;
    Real k=sqrt2*I;
    // TODO: replace Sk calculation with function pointer for efficiency?
    unsigned Sk= spectrum == DISCRETE ?  kval[I][I] : (unsigned)(k-0.5);
    Real Wall=abs2(wi[I])+abs2(wim[I]);
    S[Sk] += Complex(Wall/(k0*k),nuk(k2)*Wall);
    
    const unsigned stop=I;
    for(unsigned j=1; j < stop; ++j) {
      k2=(I2+j*j);
      k=sqrt((Real) k2);
      Sk= spectrum == DISCRETE ? kval[I][j] : (unsigned)(k-0.5);
      Wall=abs2(wi[j])+abs2(wim[j])+abs2(w[xorigin-j][I])+abs2(w[xorigin+j][I]);
      S[Sk] += Complex(Wall/(k0*k),nuk(I2+j*j)*Wall);
    }
  }

  // xorigin case
  vector wi=w[xorigin];
  for(unsigned j=1; j < my; ++j) {
    unsigned Sk= spectrum == DISCRETE ? kval[j][0] : j-1;
    Real w2=abs2(wi[j]);
    S[Sk] += Complex(w2/(k0*j),w2*nuk(j*j));
  }
  
  // bottom right
  for(unsigned i=xorigin+1; i < Nx; ++i) {
    unsigned I=i-xorigin;
    unsigned Sk= spectrum == DISCRETE ? kval[I][0] : I-1;
    Real w2=abs2(w[i][0]);
    S[Sk] += Complex(w2/(k0*I),w2*nuk(I*I));
  }
}

void DNSBase::Stochastic(const vector2&Y, double, double dt)
{
  w.Set(Y[OMEGA]);
  Set(T,Y[TRANSFER]);
  Forcing->Force(w,T,dt);
}

//***** DNSBase Output routines *****//

void DNSBase::Initialize()
{
  fevt << "# t\tE\tZ\tP" << endl;
  setcount();
}

void DNSBase::setcount()
{
  switch(spectrum) {
  case NOSPECTRUM:
    break;
  case UNINTERP: // uninterpolated spectrum
    for(unsigned i=0; i < nshells; i++)
      count[i]=0;
    for(unsigned i=0; i < Nx; i++) {
      int I=(int) i-(int) xorigin;
      int I2=I*I;
      for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
	count[(unsigned)(sqrt(I2+j*j)-0.5)]++;
      }
    }
    break;
 case INTERP: // interpolated spectrum
   msg(ERROR,"interpolated spectrum not yet implemented");
   // FIXME
   break;
 case DISCRETE: // discrete spectrum
   for(unsigned i=0; i < nshells; i++)
     count[i]=0;
   for(unsigned i=0; i < Nx; i++) {
     unsigned I= xorigin > i ? xorigin-i : i-xorigin;
     unsigned I2=I*I;
     for(unsigned j= i < xorigin?  1 : 0; j < my; ++j) {
       unsigned r2=I2+j*j;
       for(unsigned k=0; k < R2.Size(); ++k) {
	 if(r2 == R2[k]) {
	   count[k]++;
	   break;
	 }
       }
     }
   }
   break;
  }
}

void DNSBase::OutFrame(int)
{
  //  w.Set(Y[OMEGA]); // FIXME
  unsigned int Nx0=Nx+xpad;
  unsigned int Ny0=Ny+ypad;
  unsigned int offset=Nx0/2-mx+1;
  for(unsigned int i=0; i < Nx; ++i) {
    unsigned int I=i+offset;
    for(unsigned int j=0; j < my; j++)
      buffer(I,j)=w(i,j);
  }

  Padded->pad(buffer);
  Padded->backwards(buffer,true);

  fw << 1 << Ny0 << Nx0;
  for(int j=Ny0-1; j >= 0; j--) {
    for(unsigned i=0; i < Nx0; i++) {
      fw << (float) wr(i,j);
    }
  }
  fw.flush();
}

void DNSBase::ComputeInvariants(const array2<Complex> &w, 
				Real& E, Real& Z, Real& P)
{
  // TODO: optimize?
  E=Z=P=0.0;
  for(unsigned i=0; i < Nx; i++) {
    int I=(int) i-(int) xorigin;
    int I2=I*I;
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      Real w2=abs2(wi[j]);
      Z += w2;
      Real k2=(I2+j*j)*k02;
      E += w2/k2;
      P += k2*w2;
    }
  }
}

void DNSBase::FinalOutput()
{
  Real E,Z,P;
  ComputeInvariants(w,E,Z,P);
  cout << endl;
  cout << "Energy = " << E << newl;
  cout << "Enstrophy = " << Z << newl;
  cout << "Palenstrophy = " << P << newl;
}

static inline Real H(Real x)
{
  return (x > 0.0 ? x : 0.0);
}
 
static inline Real J(Real x, Real y1, Real rinv, Real r2)
{
  Real X=x*rinv;
  if(X > 1.0) X=1.0;
  return 0.5*r2*(asin(X)+X*sqrt(1.0-X*X))-x*y1;
}

static inline Real ShellArea1(Real x1, Real x2, Real y1, Real y2, Real r)
{
  Real r2=r*r;
  if(x1*x1+y1*y1 >= r2) return 0.0;
  if(x1 >= r) return 0.0;
  if(x2*x2+y2*y2 <= r2) return (x2-x1)*(y2-y1);
  
  Real rinv=1.0/r;
  Real ry1=sqrt(r2-y1*y1);
  Real ry2=sqrt(H(r2-y2*y2));
  return H(min(x2,ry2)-x1)*(y2-y1)+
    J(min(x2,ry1),y1,rinv,r2)-J(max(x1,ry2),y1,rinv,r2);
}

// Restricted to first quadrant
static inline Real ShellArea1(Real x1, Real x2, Real y1, Real y2, Real r1,
			      Real r2) 
{
  return ShellArea1(x1,x2,y1,y2,r2)-ShellArea1(x1,x2,y1,y2,r1);
}

Real ShellArea(Real x1, Real x2, Real y1, Real y2, Real r1, Real r2)
{
  if(x1 <= 0 && x2 <= 0) {x1=-x1; x2=-x2;}
  if(y1 <= 0 && y2 <= 0) {y1=-y1; y2=-y2;}
  
  if(x1 > x2) swap(x1,x2);
  if(y1 > y2) swap(y1,y2);

  if(x1 >= 0 && y1 >= 0) return ShellArea1(x1,x2,y1,y2,r1,r2);
  if(x1 < 0 && y1 < 0) {
    return ShellArea1(0,-x1,0,-y1,r1,r2)+ShellArea1(0,-x1,0,y2,r1,r2)
      + ShellArea1(0,x2,0,-y1,r1,r2)+ShellArea1(0,x2,0,y2,r1,r2);
  }
  if(x1 < 0)
    return ShellArea1(0,-x1,y1,y2,r1,r2)+ShellArea1(0,x2,y1,y2,r1,r2);
  return ShellArea1(x1,x2,0,-y1,r1,r2)+ShellArea1(x1,x2,0,y2,r1,r2);
}
