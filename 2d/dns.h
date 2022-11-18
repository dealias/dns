#ifndef __dnsbase_h__
#define __dnsbase_h__ 1

#include "options.h"
#include "kernel.h"
#include "Array.h"
#include "array2h.h"
#include "fftw++.h"
#include "convolve.h"
#include "Forcing.h"
#include "InitialCondition.h"
#include "Conservative.h"
#include "Exponential.h"
#include <sys/stat.h> // On Sun computers this must come after xstream.h

using namespace Array;
using namespace fftwpp;

using std::ostringstream;
using fftwpp::twopi;

typedef array1<Var>::opt vector;
typedef array1<Real>::opt rvector;

extern unsigned spectrum;

extern int pH;
extern int pL;

class DNSBase {
protected:
  // Vocabulary:
  unsigned Nx;
  unsigned Ny;
  Real nuH,nuL;
  Real kH2,kL2;

  // Contiguous: TRANSFERE,TRANSFERZ,EPS,ETA,ZETA,DISSIPATIONE,DISSIPATIONZ
  //
  enum Field {OMEGA,TRANSFERE,TRANSFERZ,EPS,ETA,ZETA,DISSIPATIONE,
              DISSIPATIONZ,EK};

  int mx,my; // size of data arrays
  int my1;   // x stride

  array2h<Complex> w; // Vorticity field
  array2<Real> wr; // Inverse Fourier transform of vorticity field

  int tcount;
  unsigned fcount;

  unsigned nmode;
  unsigned nshells;  // Number of spectral shells

  Array2<Complex> u,v,V;
  array2h<Complex> S;
  array2<Complex> buffer;
  Complex *F[2];
  Complex **block;

  ConvolutionHermitian2 *Convolution;
  crfft2d *Backward;

  ifstream ftin;
  oxstream fek,fw,fekvk,ftransfer;
  ofstream ft,fevt;

  uvector count;

  typedef array1<Var>::opt vector;

  vector TE,TZ; // Energy and enstrophy transfers
  vector Eps,Eta,Zeta; // Energy, enstrophy, and palenstrophy injection rates
  vector DE,DZ; // Energy and enstrophy dissipation rates
  vector E; // Energy spectrum

  Real Energy,Enstrophy,Palinstrophy;
  array2h<Real> k2inv;

public:
  array1<Real> nu;

  void Initialize() {
    fevt << "# t\tE\tZ\tP" << endl;
  }

  void InitialConditions() {
    w[0][0]=0.0; // Enforce no mean flow
    Loop(Initw(this),InitializeValue(this));
  }

  void SetParameters() {
    setcount();
    fcount=0;

    Forcing->Init();

    Loop(InitNone(this),ForcingCount(this));

    fcount *= 2; // Account for Hermitian conjugate modes.

    Forcing->Init(fcount);

    k2inv.Allocate(mx,my);
    Loop(InitNone(this),K2inv(this,sqrt(Convolution->scale)));
  }

  virtual void setcount() {
#pragma omp parallel for num_threads(threads)
    for(unsigned i=0; i < nshells; i++)
      count[i]=0;

    if(spectrum)
      Loop(InitNone(this),Count(this));
  }

  virtual void OutEnergies() {
    fek << mx << my;
    Loop(Initw(this),OutEk(this));
  }

  void FinalOutput() {
    Real E,Z,P;
    ComputeInvariants(w,E,Z,P);
    cout << endl;
    cout << "Energy = " << E << newl;
    cout << "Enstrophy = " << Z << newl;
    cout << "Palinstrophy = " << P << newl;
  }

  void OutFrame(int it) {
    Loop(Initw(this),wtoV(this));
    V(0,0)=0.0;

    fftwpp::HermitianSymmetrizeX(mx,my,mx,V);

// Zero Nyquist modes.
    for(int j=0; j < my; ++j)
      V(j)=0.0;

    Backward->fft0(V);

    fw << 1 << 2*my << Nx+1;
    for(int j=2*my-1; j >= 0; j--)
      for(unsigned i=0; i <= Nx; i++)
        fw << (float) wr(i,j);
    fw.flush();
  }

  class FETL {
    DNSBase *b;
    const vector& TE,TZ,Eps,Eta,Zeta,DE,DZ,E;
    unsigned int l;

  public:
    FETL(DNSBase *b) : b(b), TE(b->TE), TZ(b->TZ),
                       Eps(b->Eps), Eta(b->Eta), Zeta(b->Zeta),
                       DE(b->DE), DZ(b->DZ), E(b->E), l(0) {}
    inline void operator()(const vector& wi, const vector& Si, int i, int j) {
      unsigned k2=i*i+j*j;
      Real k=sqrt(k2);
      unsigned index=(unsigned)(k-0.5);
      Complex wij=wi[j];
      Real w2=abs2(wij);
      Complex& Sij=Si[j];
      Real transfer=realproduct(Sij,wij);
      Real eta=Forcing->Force(wij,Sij,i,j);
      Real kinv=1.0/k;
      Real kinv2=kinv*kinv;
      Nu nuk2=b->nu[l++];
      Real nuk2Z=nuk2*w2;
      TE[index] += kinv2*transfer;
      TZ[index] += transfer;
      Eps[index] += kinv2*eta;
      Eta[index] += eta;
      Zeta[index] += k2*eta;
      DE[index] += kinv2*nuk2Z;
      DZ[index] += nuk2Z;
      E[index] += kinv*w2;
      Sij -= nuk2*wij;
    }
  };

  class FTL {
    DNSBase *b;
    const vector& TE,TZ,Eps,Eta,Zeta,DE,DZ;
    unsigned int l;

  public:
    FTL(DNSBase *b) : b(b), TE(b->TE), TZ(b->TZ),
                      Eps(b->Eps), Eta(b->Eta), Zeta(b->Zeta),
                      DE(b->DE), DZ(b->DZ), l(0) {}
    inline void operator()(const vector& wi, const vector& Si, int i, int j) {
      unsigned k2=i*i+j*j;
      Real k=sqrt(k2);
      unsigned index=(unsigned)(k-0.5);
      Complex wij=wi[j];
      Complex& Sij=Si[j];
      Real transfer=realproduct(Sij,wij);
      Real eta=Forcing->Force(wij,Sij,i,j);
      Real kinv2=1.0/k2;
      Nu nuk2=b->nu[l++];
      Real nuk2Z=nuk2*abs2(wij);
      TE[index] += kinv2*transfer;
      TZ[index] += transfer;
      Eps[index] += kinv2*eta;
      Eta[index] += eta;
      Zeta[index] += k2*eta;
      DE[index] += kinv2*nuk2Z;
      DZ[index] += nuk2Z;
      Sij -= nuk2*wij;
    }
  };

  class FET {
    DNSBase *b;
    const vector& TE,TZ,Eps,Eta,Zeta,DE,DZ,E;
    unsigned int l;

  public:
    FET(DNSBase *b) : b(b), TE(b->TE), TZ(b->TZ),
                      Eps(b->Eps), Eta(b->Eta), Zeta(b->Zeta),
                      DE(b->DE), DZ(b->DZ), E(b->E), l(0) {}
    inline void operator()(const vector& wi, const vector& Si, int i, int j) {
      unsigned k2=i*i+j*j;
      Real k=sqrt(k2);
      unsigned index=(unsigned)(k-0.5);
      Complex wij=wi[j];
      Real w2=abs2(wij);
      Complex& Sij=Si[j];
      Real transfer=realproduct(Sij,wij);
      Real eta=Forcing->Force(wij,Sij,i,j);
      Real kinv=1.0/k;
      Real kinv2=kinv*kinv;
      Real nuk2Z=b->nu[l++]*w2;
      TE[index] += kinv2*transfer;
      TZ[index] += transfer;
      Eps[index] += kinv2*eta;
      Eta[index] += eta;
      Zeta[index] += k2*eta;
      DE[index] += kinv2*nuk2Z;
      DZ[index] += nuk2Z;
      E[index] += kinv*w2;
    }
  };

  class FE {
    DNSBase *b;
    const vector& E;

  public:
    FE(DNSBase *b) : b(b), E(b->E) {}
    inline void operator()(const vector& wi, const vector& Si, int i, int j) {
      unsigned k2=i*i+j*j;
      Real k=sqrt(k2);
      unsigned index=(unsigned)(k-0.5);
      E[index] += abs2(wi[j])/k;
    }
  };

  class FL {
    DNSBase *b;
    unsigned int l;

  public:
    FL(DNSBase *b) : b(b), l(0) {}
    inline void operator()(const vector& wi, const vector& Si, int i, int j) {
      Complex wij=wi[j];
      Forcing->Force(wij,Si[j],i,j);
      Si[j] -= b->nu[l++]*wij;
    }
  };

  class ForceStochastic {
    const vector& Eps,Eta,Zeta;
  public:
    ForceStochastic(DNSBase *b) : Eps(b->Eps), Eta(b->Eta), Zeta(b->Zeta) {}
    inline void operator()(const vector& wi, const vector&, int i, int j) {
      unsigned k2=i*i+j*j;
      Real k=sqrt(k2);
      unsigned index=(unsigned)(k-0.5);
      double eta=Forcing->ForceStochastic(wi[j],i,j);
      Eps[index] += eta/k2;
      Eta[index] += eta;
      Zeta[index] += k2*eta;
    }
  };

  class ForceStochasticNO {
  public:
    ForceStochasticNO(DNSBase *b) {}
    inline void operator()(const vector& wi, const vector&, int i, int j) {
      Forcing->ForceStochastic(wi[j],i,j);
    }
  };

  class InitializeValue {
  public:
    InitializeValue(DNSBase *b) {}
    inline void operator()(const vector& wi, const vector& Si, int i, int j) {
      wi[j]=InitialCondition->Value(i,j);
    }
  };

  class ForcingCount {
    DNSBase *b;
  public:
    ForcingCount(DNSBase *b) : b(b) {}
    inline void operator()(const vector&, const vector&, int i, int j) {
      if(Forcing->active(i,j)) {
        ++b->fcount;
      }
    }
  };

  class K2inv {
    DNSBase *b;
    double scale;
  public:
    K2inv(DNSBase *b, double scale) : b(b), scale(scale) {}
    inline void operator()(const vector&, const vector&, int i, int j) {
      b->k2inv(i,j)=scale/(i*i+j*j);
    }
  };

  class wtoV {
    DNSBase *b;
  public:
    wtoV(DNSBase *b) : b(b) {}
    inline void operator()(const vector& wi, const vector&, int i, int j) {
      b->V(i,j)=wi[j];
    }
  };

  class OutEk {
    DNSBase *b;
  public:
    OutEk(DNSBase *b) : b(b) {}
    inline void operator()(const vector& wi, const vector&, int i, int j) {
      b->fek << 0.5*abs2(wi[j])*b->k2inv(i,j);
    }
  };

  class Invariants {
    DNSBase *b;
    Real &Energy,&Enstrophy,&Palinstrophy;
  public:
    Invariants(DNSBase *b) : b(b), Energy(b->Energy), Enstrophy(b->Enstrophy),
                             Palinstrophy(b->Palinstrophy) {}
    inline void operator()(const vector& wi, const vector& Si, int i, int j) {
      Real w2=abs2(wi[j]);
      Enstrophy += w2;
      unsigned k2=i*i+j*j;
      Energy += w2/k2;
      Palinstrophy += k2*w2;
    }
  };

  class InitwS {
    DNSBase *b;
  public:
    InitwS(DNSBase *b) : b(b) {}
    inline void operator()(vector& wi, vector& Si, int i) {
      Dimension(wi,b->w[i]);
      Dimension(Si,b->S[i]);
    }
  };

  class Initw {
    DNSBase *b;
  public:
    Initw(DNSBase *b) : b(b) {}
    inline void operator()(vector& wi, vector& Si, int i) {
      Dimension(wi,b->w[i]);
    }
  };

  class InitNone {
  public:
    InitNone(DNSBase *b) {}
    inline void operator()(vector& wi, vector& Si, int i) {
    }
  };

  class Count {
    DNSBase *b;
    const uvector& count;
  public:
    Count(DNSBase *b) : b(b), count(b->count) {}
    inline void operator()(const vector& wi, const vector& Si, int i, int j) {
      unsigned k2=i*i+j*j;
      Real k=sqrt(k2);
      unsigned index=(unsigned)(k-0.5);
      ++count[index];
    }
  };

  void NonLinearSource(const vector2& Src, const vector2& Y, double t) {
    w.Set(Y[OMEGA]);

    // This 2D version of the scheme of Basdevant, J. Comp. Phys, 50, 1983
    // requires only 4 FFTs per stage.
    PARALLEL(
      for(int i=-mx+1; i < 0; ++i) {
        vector wi=w[i];
        vector ui=u[i];
        vector vi=v[i];
        rvector k2invi=k2inv[i];
        for(int j=1; j < my; ++j) {
          Complex wij=wi[j];
          Real k2invij=k2invi[j];
          Real jk2inv=j*k2invij;
          Real ik2inv=i*k2invij;
          ui[j]=Complex(-wij.im*jk2inv,wij.re*jk2inv); // u
          vi[j]=Complex(wij.im*ik2inv,-wij.re*ik2inv); // v
        }
      });


    Complex *v0=v[0];
    PARALLEL(
    for(int j=0; j < my; ++j)
      v0[j]=0.0;
      );

    PARALLEL(
    for(int i=-mx+1; i < mx; ++i)
      u[i][0]=0.0;
      );

    vector wi=w[0];
    vector ui=u[0];
    rvector k2invi=k2inv[0];
    PARALLEL(
      for(int j=1; j < my; ++j) {
        Complex wij=wi[j];
        Real jk2inv=j*k2invi[j];
        ui[j]=Complex(-wij.im*jk2inv,wij.re*jk2inv);
      });

    PARALLEL(
      for(int i=1; i < mx; ++i) {
        vector wi=w[i];
        vector ui=u[i];
        vector vi=v[i];
        rvector k2invi=k2inv[i];
        Complex wij=wi[0];
        Real ik2inv=i*k2invi[0];
        Complex V=Complex(wij.im*ik2inv,-wij.re*ik2inv);
        vi[0]=V;
        v[-i][0]=conj(V);
        for(int j=1; j < my; ++j) {
          Complex wij=wi[j];
          Real k2invij=k2invi[j];
          Real jk2inv=j*k2invij;
          Real ik2inv=i*k2invij;
          ui[j]=Complex(-wij.im*jk2inv,wij.re*jk2inv);
          vi[j]=Complex(wij.im*ik2inv,-wij.re*ik2inv);
        }
      })


     // Zero Nyquist modes.
    for(int j=0; j < my; ++j)
      u[-mx][j]=0.0;
    for(int j=0; j < my; ++j)
      v[-mx][j]=0.0;

    Convolution->convolveRaw(F);

    S.Set(Src[OMEGA]);

    PARALLEL(
    for(int i=-mx+1; i < mx; ++i) {
      Real i2=i*i;
      vector ui=u[i];
      vector vi=v[i];
      vector Si=S[i];
      for(int j=i <= 0; j < my; ++j) {
        Si[j]=i*j*ui[j]+(i2-j*j)*vi[j];
      }
    });


#if 0
    Real sum=0.0;
    for(int i=-mx+1; i < mx; ++i) {
      vector wi=w[i];
      for(int j=i <= 0; j < my; ++j) {
        Complex wij=wi[j];
//        sum += (S[i][j]*conj(wij)).re;
        sum += (S[i][j]*conj(wij)/(i*i+j*j)).re;
      }
    }
    cout << sum << endl;
    cout << endl;
#endif
  }

  void Init(vector& T, const vector& Src) {
    Set(T,Src);
#pragma omp parallel for num_threads(threads)
    for(unsigned K=0; K < nshells; K++)
      T[K]=0.0;
  }

  void ConservativeSource(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
    if(spectrum) {
      Init(TE,Src[TRANSFERE]);
      Init(TZ,Src[TRANSFERZ]);
      Init(Eps,Src[EPS]);
      Init(Eta,Src[ETA]);
      Init(Zeta,Src[ZETA]);
      Init(DE,Src[DISSIPATIONE]);
      Init(DZ,Src[DISSIPATIONZ]);
      Compute(FTL(this),Src,Y);
    }
    else
      Compute(FL(this),Src,Y);
  }

  void NonConservativeSource(const vector2& Src, const vector2& Y, double t) {
    if(spectrum) {
      Init(E,Src[EK]);
      Compute(FE(this),Src,Y);
    }
  }

  void ExponentialSource(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
    if(spectrum) {
      Init(TE,Src[TRANSFERE]);
      Init(TZ,Src[TRANSFERZ]);
      Init(Eps,Src[EPS]);
      Init(Eta,Src[ETA]);
      Init(Zeta,Src[ZETA]);
      Init(DE,Src[DISSIPATIONE]);
      Init(DZ,Src[DISSIPATIONZ]);
      Init(E,Src[EK]);
      Compute(FET(this),Src,Y);
    }
  }

  void Source(const vector2& Src, const vector2& Y, double t) {
    NonLinearSource(Src,Y,t);
    if(spectrum) {
      Init(TE,Src[TRANSFERE]);
      Init(TZ,Src[TRANSFERZ]);
      Init(Eps,Src[EPS]);
      Init(Eta,Src[ETA]);
      Init(Zeta,Src[ZETA]);
      Init(DE,Src[DISSIPATIONE]);
      Init(DZ,Src[DISSIPATIONZ]);
      Init(E,Src[EK]);
      Compute(FETL(this),Src,Y);
    } else
      Compute(FL(this),Src,Y);
  }

  template<class S, class T>
  void Loop(S init, T fcn)
  {
    vector wi,Si;
    for(int i=-mx+1; i < mx; ++i) {
      init(wi,Si,i);
      for(int j=i <= 0; j < my; ++j) { // start with j=1 if i <= 0
        fcn(wi,Si,i,j);
      }
    }
  }

  template<class T>
  void Compute(T fcn, const vector2& Src, const vector2& Y)
  {
    S.Set(Src[OMEGA]);
    w.Set(Y[OMEGA]);

    Loop(InitwS(this),fcn);
  }

  void Stochastic(const vector2&Y, double, double dt)
  {
    if(!Forcing->Stochastic(dt)) return;
    w.Set(Y[OMEGA]);

    if(spectrum == 0) {
      Loop(Initw(this),ForceStochasticNO(this));
    } else {
      Set(Eps,Y[EPS]);
      Set(Eta,Y[ETA]);
      Set(Zeta,Y[ZETA]);
      Loop(Initw(this),ForceStochastic(this));
    }
  }

  Nu LinearCoeff(unsigned l) {
    return nu[l];
  }

  Real nuk(double k2) {
    Real diss=0.0;
    if(k2 < kL2) diss += nuL*pow(k2,pL);
    if(k2 >= kH2) diss += nuH*pow(k2,pH);
    return diss;
  }

  virtual void ComputeInvariants(const array2h<Complex> &w, Real& E, Real& Z,
                                 Real& P) {
    Energy=Enstrophy=Palinstrophy=0.0;

    Loop(Initw(this),Invariants(this));

    E=Energy;
    Z=Enstrophy;
    P=Palinstrophy;
  }

  virtual Real getSpectrum(unsigned i) {
    double c=count[i];
    return c > 0 ? E[i].re*twopi/c : 0.0;
  }
  Real TE_(unsigned i) {return TE[i].re;}
  Real TZ_(unsigned i) {return TZ[i].re;}
  Real Eps_(unsigned i) {return Eps[i].re;}
  Real Eta_(unsigned i) {return Eta[i].re;}
  Real Zeta_(unsigned i) {return Zeta[i].re;}
  Real DE_(unsigned i) {return DE[i].re;}
  Real DZ_(unsigned i) {return DZ[i].re;}

  Real kb(unsigned i) {return i+0.5;}
  Real kc(unsigned i) {return i+1;}
};

extern InitialConditionBase *InitialCondition;

#endif
