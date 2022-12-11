#ifndef __dnsbase_h__
#define __dnsbase_h__ 1

#include "options.h"

typedef ptrdiff_t Int;
typedef size_t uInt;

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
#include "seconds.h"

using namespace Array;
using namespace fftwpp;

using std::ostringstream;
using fftwpp::twopi;

typedef array1<Var>::opt vector;
typedef array1<Real>::opt rvector;
typedef array1<Nu>::opt nuvector;

extern double sum;

extern uInt spectrum;

extern int pH;
extern int pL;

class DNSBase {
protected:
  // Vocabulary:
  uInt Nx;
  uInt Ny;
  Real nuH,nuL;
  Real kH2,kL2;

  // Contiguous: TRANSFERE,TRANSFERZ,EPS,ETA,ZETA,DISSIPATIONE,DISSIPATIONZ
  //
  enum Field {OMEGA,TRANSFERE,TRANSFERZ,EPS,ETA,ZETA,DISSIPATIONE,
              DISSIPATIONZ,EK};

  Int mx,my; // size of data arrays
  Int my1;   // x stride

  array2h<Complex> w; // Vorticity field
  array2<Real> wr; // Inverse Fourier transform of vorticity field

  uInt tcount;
  uInt fcount;

  uInt nmode;
  uInt nshells;  // Number of spectral shells

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

  array2h<Real> k2inv;
  array2h<Real> nu; // Linear dissipation

  Real Energy,Enstrophy,Palinstrophy;

public:

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

    nu.Allocate(mx,my);
    Loop(Initnu(this),Linearity(this));
  }

  virtual void setcount() {
    PARALLELIF(
      nshells > threshold,
      for(uInt i=0; i < nshells; i++)
        count[i]=0;
      );

    if(spectrum)
      Loop(InitNone(this),Count(this));
  }

  virtual void OutEnergies() {
    fek << mx << my;
    Loop(Initw(this),OutEk(this));
  }

  void FinalOutput();

  void OutFrame(uInt it) {
    Loop(Initw(this),wtoV(this));
    V(0,0)=0.0;

    fftwpp::HermitianSymmetrizeX(mx,my,mx,V);

// Zero Nyquist modes.
    PARALLELIF(
      my > (Int) threshold,
      for(Int j=0; j < my; ++j)
        V(j)=0.0;
      );

    Backward->fft0(V);

    fw << 1 << 2*my << Nx+1;
    for(Int j=2*my-1; j >= 0; j--)
      for(uInt i=0; i <= Nx; i++)
        fw << (float) wr(i,j);
    fw.flush();
  }

  class FETL {
    DNSBase *b;
    const vector& TE,TZ,Eps,Eta,Zeta,DE,DZ,E;

  public:
    FETL(DNSBase *b) : b(b), TE(b->TE), TZ(b->TZ),
                       Eps(b->Eps), Eta(b->Eta), Zeta(b->Zeta),
                       DE(b->DE), DZ(b->DZ), E(b->E) {}
    inline void operator()(const vector& wi, const vector& Si,
                           const nuvector &nui, Int i, Int j) {
      uInt k2=i*i+j*j;
      Real k=sqrt(k2);
      uInt index=(uInt)(k-0.5);
      Complex wij=wi[j];
      Real w2=abs2(wij);
      Complex& Sij=Si[j];
      Real transfer=realproduct(Sij,wij);
      Real eta=Forcing->Force(wij,Sij,i,j);
      Real kinv=1.0/k;
      Real kinv2=kinv*kinv;
      Nu nuk2=nui[j];
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

  public:
    FTL(DNSBase *b) : b(b), TE(b->TE), TZ(b->TZ),
                      Eps(b->Eps), Eta(b->Eta), Zeta(b->Zeta),
                      DE(b->DE), DZ(b->DZ) {}
    inline void operator()(const vector& wi, const vector& Si,
                           const nuvector &nui, Int i, Int j) {
      uInt k2=i*i+j*j;
      Real k=sqrt(k2);
      uInt index=(uInt)(k-0.5);
      Complex wij=wi[j];
      Complex& Sij=Si[j];
      Real transfer=realproduct(Sij,wij);
      Real eta=Forcing->Force(wij,Sij,i,j);
      Real kinv2=1.0/k2;
      Nu nuk2=nui[j];
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

  public:
    FET(DNSBase *b) : b(b), TE(b->TE), TZ(b->TZ),
                      Eps(b->Eps), Eta(b->Eta), Zeta(b->Zeta),
                      DE(b->DE), DZ(b->DZ), E(b->E) {}
    inline void operator()(const vector& wi, const vector& Si,
                           const nuvector &nui, Int i, Int j) {
      uInt k2=i*i+j*j;
      Real k=sqrt(k2);
      uInt index=(uInt)(k-0.5);
      Complex wij=wi[j];
      Real w2=abs2(wij);
      Complex& Sij=Si[j];
      Real transfer=realproduct(Sij,wij);
      Real eta=Forcing->Force(wij,Sij,i,j);
      Real kinv=1.0/k;
      Real kinv2=kinv*kinv;
      Real nuk2Z=nui[j]*w2;
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
    inline void operator()(const vector& wi, const vector&,
                           const nuvector &, Int i, Int j) {
      uInt k2=i*i+j*j;
      Real k=sqrt(k2);
      uInt index=(uInt)(k-0.5);
      E[index] += abs2(wi[j])/k;
    }
  };

  class FL {
    DNSBase *b;

  public:
    FL(DNSBase *b) : b(b) {}
    inline void operator()(const vector& wi, const vector& Si,
                           const nuvector &nui, Int i, Int j) {
      Complex wij=wi[j];
      Forcing->Force(wij,Si[j],i,j);
      Si[j] -= nui[j]*wij;
    }
  };

  class ForceStochastic {
    const vector& Eps,Eta,Zeta;
  public:
    ForceStochastic(DNSBase *b) : Eps(b->Eps), Eta(b->Eta), Zeta(b->Zeta) {}
    inline void operator()(const vector& wi, const vector& Si,
                           const nuvector &, Int i, Int j) {
      uInt k2=i*i+j*j;
      Real k=sqrt(k2);
      uInt index=(uInt)(k-0.5);
      double eta=Forcing->ForceStochastic(wi[j],i,j);
      Eps[index] += eta/k2;
      Eta[index] += eta;
      Zeta[index] += k2*eta;
    }
  };

  class ForceStochasticNO {
  public:
    ForceStochasticNO(DNSBase *b) {}
    inline void operator()(const vector& wi, const vector&,
                           const nuvector &, Int i, Int j) {
      Forcing->ForceStochastic(wi[j],i,j);
    }
  };

  class InitializeValue {
  public:
    InitializeValue(DNSBase *b) {}
    inline void operator()(const vector& wi, const vector&,
                           const nuvector &, Int i, Int j) {
      wi[j]=InitialCondition->Value(i,j);
    }
  };

  class ForcingCount {
    DNSBase *b;
  public:
    ForcingCount(DNSBase *b) : b(b) {}
    inline void operator()(const vector&, const vector&,
                           const nuvector &, Int i, Int j) {
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
    inline void operator()(const vector&, const vector&,
                           const nuvector &, Int i, Int j) {
      b->k2inv(i,j)=scale/(i*i+j*j);
    }
  };

  class Linearity {
    DNSBase *b;
  public:
    Linearity(DNSBase *b) : b(b) {}
    inline void operator()(const vector&, const vector&,
                           const nuvector &nui, Int i, Int j) {
      nui[j]=b->nuk(i*i+j*j);
    }
  };

  class wtoV {
    DNSBase *b;
  public:
    wtoV(DNSBase *b) : b(b) {}
    inline void operator()(const vector&wi, const vector&,
                           const nuvector &nui, Int i, Int j) {
      b->V(i,j)=wi[j];
    }
  };

  class OutEk {
    DNSBase *b;
  public:
    OutEk(DNSBase *b) : b(b) {}
    inline void operator()(const vector&wi, const vector&,
                           const nuvector &, Int i, Int j) {
      b->fek << 0.5*abs2(wi[j])*b->k2inv(i,j);
    }
  };

  class Invariants {
    Real Energy,Enstrophy,Palinstrophy;
  public:
    Invariants() : Energy(0.0), Enstrophy(0.0), Palinstrophy(0.0) {}
    inline void operator()(const vector& wi, Int i, Int j) {
      Real w2=abs2(wi[j]);
      Enstrophy += w2;
      uInt k2=i*i+j*j;
      Energy += w2/k2;
      Palinstrophy += k2*w2;
    }

    void reduce(Real& E, Real& Z, Real& P) {
      E += Energy;
      Z += Enstrophy;
      P += Palinstrophy;
    }
  };

  class InitAll {
    DNSBase *b;
  public:
    InitAll(DNSBase *b) : b(b) {}
    inline void operator()(vector& wi, vector& Si, nuvector &nui, Int i) {
      Dimension(wi,b->w[i]);
      Dimension(Si,b->S[i]);
      Dimension(nui,b->nu[i]);
    }
  };

  class Initw {
    DNSBase *b;
  public:
    Initw(DNSBase *b) : b(b) {}
    inline void operator()(vector& wi, vector&, nuvector&, Int i) {
      Dimension(wi,b->w[i]);
    }
  };

  class Initnu {
    DNSBase *b;
  public:
    Initnu(DNSBase *b) : b(b) {}
    inline void operator()(vector&, vector&, nuvector& nui, Int i) {
      Dimension(nui,b->nu[i]);
    }
  };

  class InitNone {
  public:
    InitNone(DNSBase *b) {}
    inline void operator()(vector&, vector&, nuvector&, Int) {
    }
  };

  class Count {
    DNSBase *b;
    const uvector& count;
  public:
    Count(DNSBase *b) : b(b), count(b->count) {}
    inline void operator()(const vector&, const vector&,
                           const nuvector &, Int i, Int j) {
      uInt k2=i*i+j*j;
      Real k=sqrt(k2);
      uInt index=(uInt)(k-0.5);
      ++count[index];
    }
  };

  void NonLinearSource(const vector2& Src, const vector2& Y, double t) {
    w.Set(Y[OMEGA]);

    // This 2D version of the scheme of Basdevant, J. Comp. Phys, 50, 1983
    // requires only 4 FFTs per stage.
    PARALLELIF(
      mx*my > (Int) threshold,
      for(Int i=-mx+1; i < 0; ++i) {
        vector wi=w[i];
        vector ui=u[i];
        vector vi=v[i];
        rvector k2invi=k2inv[i];
        for(Int j=1; j < my; ++j) {
          Complex wij=wi[j];
          Real k2invij=k2invi[j];
          Real jk2inv=j*k2invij;
          Real ik2inv=i*k2invij;
          ui[j]=Complex(-wij.im*jk2inv,wij.re*jk2inv); // u
          vi[j]=Complex(wij.im*ik2inv,-wij.re*ik2inv); // v
        }
      });


    Complex *v0=v[0];

    PARALLELIF(
      my > (Int) threshold,
    for(Int j=0; j < my; ++j)
      v0[j]=0.0;
      );

    PARALLELIF(
      2*mx > (Int) threshold,
      for(Int i=-mx+1; i < mx; ++i)
        u[i][0]=0.0;
      );

    vector wi=w[0];
    vector ui=u[0];
    rvector k2invi=k2inv[0];
    PARALLELIF(
      my > (Int) threshold,
      for(Int j=1; j < my; ++j) {
        Complex wij=wi[j];
        Real jk2inv=j*k2invi[j];
        ui[j]=Complex(-wij.im*jk2inv,wij.re*jk2inv);
      });

    PARALLELIF(
      mx*my > (Int) threshold,
      for(Int i=1; i < mx; ++i) {
        vector wi=w[i];
        vector ui=u[i];
        vector vi=v[i];
        rvector k2invi=k2inv[i];
        Complex wij=wi[0];
        Real ik2inv=i*k2invi[0];
        Complex V=Complex(wij.im*ik2inv,-wij.re*ik2inv);
        vi[0]=V;
        v[-i][0]=conj(V);
        for(Int j=1; j < my; ++j) {
          Complex wij=wi[j];
          Real k2invij=k2invi[j];
          Real jk2inv=j*k2invij;
          Real ik2inv=i*k2invij;
          ui[j]=Complex(-wij.im*jk2inv,wij.re*jk2inv);
          vi[j]=Complex(wij.im*ik2inv,-wij.re*ik2inv);
        }
      })


     // Zero Nyquist modes.
    PARALLELIF(
      my > (Int) threshold,
      for(Int j=0; j < my; ++j)
        u[-mx][j]=0.0;
      );

    PARALLELIF(
      my > (Int) threshold,
      for(Int j=0; j < my; ++j)
        v[-mx][j]=0.0;
      );

    double t0=utils::nanoseconds();
    Convolution->convolveRaw(F);
    sum += (utils::nanoseconds()-t0)*1.0e-9;

    S.Set(Src[OMEGA]);

    PARALLELIF(
      2*mx > (Int) threshold,
      for(Int i=-mx+1; i < mx; ++i) {
        Real i2=i*i;
        vector ui=u[i];
        vector vi=v[i];
        vector Si=S[i];
        for(Int j=i <= 0; j < my; ++j) {
          Si[j]=i*j*ui[j]+(i2-j*j)*vi[j];
        }
      });

#if 0
    Real sum=0.0;
    for(Int i=-mx+1; i < mx; ++i) {
      vector wi=w[i];
      for(Int j=i <= 0; j < my; ++j) {
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
    PARALLELIF(
      nshells > threshold,
      for(uInt K=0; K < nshells; K++)
        T[K]=0.0;
      );
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
    // TODO: Allocate separate TE...E arrays for each thread

//    PARALLELIF(
//      2*mx > (Int) threshold,
      for(Int i=-mx+1; i < mx; ++i) {
        vector wi;
        vector Si;
        nuvector nui;
        init(wi,Si,nui,i);
        for(Int j=i <= 0; j < my; ++j) // start with j=1 if i <= 0
          fcn(wi,Si,nui,i,j);
      }//);
  }

  template<class T>
  void Compute(T fcn, const vector2& Src, const vector2& Y)
  {
    S.Set(Src[OMEGA]);
    w.Set(Y[OMEGA]);

    Loop(InitAll(this),fcn);
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

  Nu LinearCoeff(uInt l) {
    return nu(l);
  }

  Real nuk(double k2) {
    Real diss=0.0;
    if(k2 < kL2) diss += nuL*pow(k2,pL);
    if(k2 >= kH2) diss += nuH*pow(k2,pH);
    return diss;
  }

  virtual void ComputeInvariants(const array2h<Complex> &w) {
    Real E=0.0, Z=0.0, P=0.0;

    OMPIF(
      2*mx > (Int) threshold,
      "omp parallel for num_threads(threads) reduction(+ : E,Z,P)",
      for(Int i=-mx+1; i < mx; ++i) {
        vector wi=this->w[i];
        Invariants I;
        for(Int j=i <= 0; j < my; ++j) // start with j=1 if i <= 0
          I(wi,i,j);
        I.reduce(E,Z,P);
      });

    this->Energy=E;
    this->Enstrophy=Z;
    this->Palinstrophy=P;
  }

  virtual Real getSpectrum(uInt i) {
    double c=count[i];
    return c > 0 ? E[i].re*twopi/c : 0.0;
  }
  Real TE_(uInt i) {return TE[i].re;}
  Real TZ_(uInt i) {return TZ[i].re;}
  Real Eps_(uInt i) {return Eps[i].re;}
  Real Eta_(uInt i) {return Eta[i].re;}
  Real Zeta_(uInt i) {return Zeta[i].re;}
  Real DE_(uInt i) {return DE[i].re;}
  Real DZ_(uInt i) {return DZ[i].re;}

  Real kb(uInt i) {return i+0.5;}
  Real kc(uInt i) {return i+1;}
};

extern InitialConditionBase *InitialCondition;

#endif
