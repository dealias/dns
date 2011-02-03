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
  
  killmodes(f0);
  killmodes(w);
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

  unsigned mx=(Nx+1)/2;

  // diagonals
  unsigned stop=diagstop();
  for(unsigned I=diagstart(); I < stop; I++) {
    unsigned i=xorigin+I;
    unsigned im=xorigin-I;
    
    vector wi=w[i];
    vector wim=w[im];
    vector fi=f0[i];
    vector fim=f0[im];

    Real k=sqrt2*I;
    unsigned Sk=(this->*Sindex)(I,I,k);
    T[Sk].re += realproduct(fi[I],wi[I]) + realproduct(fim[I],wim[I]);
  }

  // main case
  for(unsigned I=1; I < mx; I++) { // FIXME: add bounds here too?
    unsigned i=xorigin+I;
    unsigned im=xorigin-I;
    unsigned I2=I*I;

    vector wi=w[i];
    vector wim=w[im];
    vector fi=f0[i];
    vector fim=f0[im];

    const unsigned stop=mainjstop(I);
    for(unsigned j=mainjstart(I); j < stop; ++j) {
      unsigned Sk=(this->*Sindex)(I,j,sqrt((I2+j*j)));
      unsigned jm=xorigin-j;
      unsigned jp=xorigin+j;
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
  stop=xoriginstop();
  for(unsigned j=xoriginstart(); j < stop; ++j) {
    unsigned Sk=(this->*Sindex)(j,0,j);
    T[Sk].re += realproduct(fi[j],wi[j]);
  }
  
  // bottom right
  stop=bottomstop();
  for(unsigned I=bottomstart(); I < stop; I++) {
    unsigned i=I+xorigin;
    unsigned Sk=(this->*Sindex)(I,0,I);
    T[Sk].re += realproduct(f0[i][0],w[i][0]);
  }
}

void DNSBase::Spectrum(vector& S, const vector& y)
{
  w.Set(y);
  for(unsigned K=0; K < nshells; K++)    
    S[K]=Complex(0.0,0.0);

  unsigned mx=(Nx+1)/2;

  // diagonals
  unsigned stop=diagstop();
  for(unsigned I=diagstart(); I < stop; I++) {
    unsigned i=xorigin+I;
    unsigned im=xorigin-I;
    unsigned I2=I*I;

    vector wi=w[i];
    vector wim=w[im];

    Real k=sqrt2*I;
    unsigned k2=2*I2;
    unsigned Sk=(this->*Sindex)(I,I,k);
    Real Wall=abs2(wi[I])+abs2(wim[I]);
    S[Sk] += Complex(Wall/(k0*k),nuk(k2)*Wall);
  }

  // main case
  for(unsigned I=1; I < mx; I++) {
    unsigned i=xorigin+I;
    unsigned im=xorigin-I;
    unsigned I2=I*I;

    vector wi=w[i];
    vector wim=w[im];

    const unsigned stop=mainjstop(I);
    for(unsigned j=mainjstart(I); j < stop; ++j) {
      unsigned k2=(I2+j*j);
      Real k=sqrt((Real) k2);
      unsigned Sk=(this->*Sindex)(I,j,k);
      Real Wall=abs2(wi[j])+abs2(wim[j])
	+abs2(w[xorigin-j][I])+abs2(w[xorigin+j][I]);
      S[Sk] += Complex(Wall/(k0*k),nuk(I2+j*j)*Wall);
    }
  }

  // xorigin case
  vector wi=w[xorigin];
  //for(unsigned j=Invis == 0 ? 1 : Invis; j < my; ++j) {
  stop=xoriginstop();
  for(unsigned j=xoriginstart(); j < stop; ++j) {
    //unsigned Sk= spectrum == RAW ? kval[j][0] : j-1;
    unsigned Sk=(this->*Sindex)(j,0,j);
    Real w2=abs2(wi[j]);
    S[Sk] += Complex(w2/(k0*j),w2*nuk(j*j));
  }
  
  // bottom right
  //for(unsigned i=Invis == 0? xorigin+1 : xorigin + Invis; i < Nx; ++i) {
  stop=bottomstop();
  for(unsigned I=bottomstart(); I < stop; I++) {
    unsigned i=I+xorigin;
    //unsigned Sk= spectrum == RAW ? kval[I][0] : I-1;
    unsigned Sk=(this->*Sindex)(I,0,I);
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

void DNSBase::setcountBINNED()
{
  for(unsigned i=0; i < Nx; i++) {
    int I=(int) i-(int) xorigin;
    int I2=I*I;
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      const Real kint=sqrt(I2+j*j);
      if(isvisible(i,j)) // FIXME
	count[(unsigned)(kint-0.5)]++;
    }
  }
}

void DNSBase::setcountRAW(unsigned lambda2)
{
  for(unsigned i=0; i < Nx; i++) {
    unsigned I= xorigin > i ? xorigin-i : i-xorigin;
    unsigned I2=I*I;
    for(unsigned j= i < xorigin?  1 : 0; j < my; ++j) {
      if(isvisible(I,j)){
	unsigned r2=lambda2*(I2+j*j);
	for(unsigned k=0; k < R2.Size(); ++k) {
	  if(r2 == R2[k]) {
	    count[k]++;
	    break;
	  }
	}
      }
    }
  }
}

void DNSBase::setcount()
{
  for(unsigned i=0; i < nshells; i++)
    count[i]=0;
  switch(spectrum) {
  case NOSPECTRUM:
    break;
  case BINNED:
    setcountBINNED();
    break;
 case INTERPOLATED:
   msg(ERROR,"interpolated spectrum not yet implemented");
   break;
 case RAW:
   setcountRAW();
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
