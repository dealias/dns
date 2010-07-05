#include "options.h"
#include "dns.h"




const double ProblemVersion=1.0;

#ifndef DEPEND
#if !COMPLEX
#error DNS requires COMPLEX=1 in options.h
#endif
#endif

const char *problem="Direct Numerical Simulation of Turbulence";

const char *method="DNS";
const char *integrator="RK5";
const char *ic="Equipartition";
//const char *linearity="Power";
const char *forcing="WhiteNoiseBanded";

// Vocabulary
Real nuH=0.0, nuL=0.0;
int pH=1;
int pL=0;
unsigned Nx=15;
unsigned Ny=15;
Real eta=0.0;
Complex force=0.0;
Real kforce=1.0;
Real deltaf=1.0;
unsigned movie=0;
unsigned rezero=0;
unsigned spectrum=1;
Real icalpha=1.0;
Real icbeta=1.0;

// other global variables:
int xpad=1;
int ypad=1;


class DNS : public DNSBase, public ProblemBase {
public:
  DNS();
  ~DNS();
  void InitialConditions();
  void Output(int);

  void IndexLimits(unsigned& start, unsigned& stop,
		   unsigned& startT, unsigned& stopT,
		   unsigned& startM, unsigned& stopM) {
    start=Start(OMEGA);
    stop=Stop(OMEGA);
    startT=Start(TRANSFER);
    stopT=Stop(TRANSFER);
    startM=Start(EK);
    stopM=Stop(EK);
  }
  void Source(const vector2& Src, const vector2& Y, double t) {
    ConservativeSource(Src,Y,t);
    NonConservativeSource(Src,Y,t);
  }

};

DNS *DNSProblem;


InitialConditionBase *InitialCondition;
ForcingBase *Forcing;

class Constant : public InitialConditionBase {
public:
  const char *Name() {return "Constant";}
  void Set(Complex *w, unsigned n) {
    for(unsigned i=0; i < n; i++)
      w[i]=Complex(icalpha,icbeta);
  }
};

class Equipartition : public InitialConditionBase {
public:
  const char *Name() {return "Equipartition";}
  void Set(Complex *w0, unsigned) {
    unsigned Nx=DNSProblem->getNx();
    unsigned my=DNSProblem->getmy();
    unsigned xorigin=DNSProblem->getxorigin();
    Real k0=DNSProblem->getk0();

    array2<Complex> w(Nx,my,w0);
    w(xorigin,0)=0;
    Real k02=k0*k0;
    for(unsigned i=0; i < Nx; i++) {
      int I=(int) i-(int) xorigin;
      int I2=I*I;
      vector wi=w[i];
      for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
	Real k2=k02*(I2+j*j);
// Distribute the enstrophy evenly between the real and imaginary components
        Real v=icalpha+icbeta*k2;
        v=v ? sqrt(0.5*k2/v) : 0.0;
	wi[j]=Complex(v,v);
      }
    }
  }
};

// forcing
class None : public ForcingBase {
};

class ConstantBanded : public ForcingBase {
public:
  const char *Name() {return "Constant Banded";}
  void Force(array2<Complex> &w, vector& T, const Complex&) {
    unsigned Nx=DNSProblem->getNx();
    unsigned my=DNSProblem->getmy();
    unsigned xorigin=DNSProblem->getxorigin();
    Real k02=DNSProblem->getk02();
    Real kmin=max(kforce-0.5*deltaf,0.0);
    Real kmin2=kmin*kmin;
    Real kmax=kforce+0.5*deltaf;
    Real kmax2=kmax*kmax;

    // TODO: only loop over modes with k in (kmin,kmax)
    for(unsigned i=0; i < Nx; i++) {
      int I=(int) i-(int) xorigin;
      int I2=I*I;
      vector wi=w[i];
      for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
	Real k2=k02*(I2+j*j);
	if(k2 > kmin2 && k2 < kmax2) {
          T[(unsigned)(sqrt(k2)-0.5)].im += realproduct(force,wi[j]);
	  wi[j] += force;
        }
      }
    }
  }
};


class WhiteNoiseBanded : public ForcingBase {
public:
  const char *Name() {return "White-Noise Banded";}
  void Force(array2<Complex> &w, const Complex& factor) {
    unsigned Nx=DNSProblem->getNx();
    unsigned my=DNSProblem->getmy();
    unsigned xorigin=DNSProblem->getxorigin();
    Real k02=DNSProblem->getk02();
    Real kmin=max(kforce-0.5*deltaf,0.0);
    Real kmin2=kmin*kmin;
    Real kmax=kforce+0.5*deltaf;
    Real kmax2=kmax*kmax;

    // TODO: only loop over modes with k in (kmin,kmax)
    Complex Factor=factor*sqrt(2.0*eta);
    for(unsigned i=0; i < Nx; i++) {
      int I=(int) i-(int) xorigin;
      int I2=I*I;
      vector wi=w[i];
      for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
	Real k2=k02*(I2+j*j);
	if(k2 > kmin2 && k2 < kmax2)
	  wi[j] += Factor;
      }
    }
  }
};

class DNSVocabulary : public VocabularyBase {
public:
  const char *Name() {return problem;}
  const char *Abbrev() {return "DNS";}
  DNSVocabulary();

  Table<InitialConditionBase> *InitialConditionTable;
  InitialConditionBase *NewInitialCondition(const char *& key) {
    return InitialConditionTable->Locate(key);
  }

  Table<ForcingBase> *ForcingTable;
  ForcingBase *NewForcing(const char *& key) {
    return ForcingTable->Locate(key);
  }
};

DNSVocabulary DNS_Vocabulary;

DNSVocabulary::DNSVocabulary()
{
  Vocabulary=this;

  VOCAB_NOLIMIT(ic,"Initial Condition");
  VOCAB(Nx,1,INT_MAX,"Number of dealiased modes in x direction");
  VOCAB(Ny,1,INT_MAX,"Number of dealiased modes in y direction");
  VOCAB(movie,0,1,"Movie flag (0=off, 1=on)");
  VOCAB(spectrum,0,1,"Spectrum flag (0=off, 1=on)");
  VOCAB(rezero,0,INT_MAX,"Rezero moments every rezero output steps for high accuracy");

  METHOD(DNS);

  InitialConditionTable=new Table<InitialConditionBase>("initial condition");
  VOCAB(icalpha,0.0,0.0,"initial condition parameter");
  VOCAB(icbeta,0.0,0.0,"initial condition parameter");
  INITIALCONDITION(Zero);
  INITIALCONDITION(Constant);
  INITIALCONDITION(Equipartition);

  VOCAB(nuH,0.0,REAL_MAX,"High-wavenumber viscosity");
  VOCAB(nuL,0.0,REAL_MAX,"Low-wavenumber viscosity");
  VOCAB(pH,0,0,"Power of Laplacian for high-wavenumber viscosity");
  VOCAB(pL,0,0,"Power of Laplacian for molecular viscosity");

  VOCAB_NOLIMIT(forcing,"Forcing type");
  ForcingTable=new Table<ForcingBase>("forcing");

  VOCAB(eta,0.0,REAL_MAX,"vorticity injection rate");
  VOCAB(force,(Complex) 0.0, (Complex) 0.0,"constant external force");
  VOCAB(kforce,0.0,REAL_MAX,"forcing wavenumber");
  VOCAB(deltaf,0.0,REAL_MAX,"forcing band width");
  FORCING(None);
  FORCING(WhiteNoiseBanded);
}

// DNSBase setup routines
DNS::DNS()
{
  DNSProblem=this;
  check_compatibility(DEBUG);
  ConservativeIntegrators(DNS_Vocabulary.IntegratorTable,this);
  ExponentialIntegrators(DNS_Vocabulary.IntegratorTable,this);
}

DNS::~DNS()
{
  fftwpp::deleteAlign(block);
}


void DNS::InitialConditions()
{
  // load vocabulary from global variables
  xpad=::xpad;
  ypad=::ypad;
  nuH=::nuH;
  nuL=::nuL;
  pH=::pH;
  pL=::pL;
  Nx=::Nx;
  Ny=::Ny;
  spectrum=::spectrum;
  if(Nx % 2 == 0 || Ny % 2 == 0) msg(ERROR,"Nx and Ny must be odd");

  k0=1.0;
  k02=k0*k0;
  mx=(Nx+1)/2;
  my=(Ny+1)/2;
  xorigin=mx-1;
  origin=xorigin*my;
  nshells=spectrum ? (unsigned) (hypot(mx-1,my-1)+0.5) : 0;

  NY[OMEGA]=Nx*my;
  NY[TRANSFER]=nshells;
  NY[EK]=nshells;

  cout << "\nGEOMETRY: (" << Nx << " X " << Ny << ")" << endl;
  cout << "\nALLOCATING FFT BUFFERS" << endl;
  size_t align=sizeof(Complex);

  Allocator(align);

  Dimension(T,nshells);

  w.Dimension(Nx,my);
  f0.Dimension(Nx,my);

  unsigned int Nxmy=Nx*my;
  unsigned int nbuf=3*Nxmy;
  unsigned int Nx0=Nx+xpad;
  unsigned int Ny0=Ny+ypad;
  int my0=Ny0/2+1;
  if(movie)
    nbuf=max(nbuf,Nx0*my0);

  block=fftwpp::ComplexAlign(nbuf);
  f1.Dimension(Nx,my,block);
  g0.Dimension(Nx,my,block+Nxmy);
  g1.Dimension(Nx,my,block+2*Nxmy);

  F[1]=f1;
  G[0]=g0;
  G[1]=g1;

  Convolution=new fftwpp::ImplicitHConvolution2(mx,my,2);

  Allocate(count,nshells);

  if(movie) {
    buffer.Dimension(Nx0,my0,block);
    wr.Dimension(Nx0,2*my0,(Real *) block);
    Padded=new fftwpp::ExplicitHConvolution2(Nx0,Ny0,mx,my,block);
  }

  InitialCondition=DNS_Vocabulary.NewInitialCondition(ic);
  w.Set(Y[OMEGA]);
  InitialCondition->Set(w,NY[OMEGA]);
  fftwpp::HermitianSymmetrizeX(mx,my,xorigin,w);

  for(unsigned i=0; i < nshells; i++)
    Y[EK][i]=0.0;

  Forcing=DNS_Vocabulary.NewForcing(forcing);

  if(dynamic && false) {
    Allocate(errmask,ny);
    for(unsigned i=0; i < ny; ++i)
      errmask[i]=1;

    array2<int> omegamask(Nx,my,(int *) errmask+(Y[OMEGA]-y));
    for(unsigned i=0; i <= xorigin; i++)
      omegamask(i,0)=0;
  }

  tcount=0;
  if(restart) {
    Real t0;
    ftin.open(Vocabulary->FileName(dirsep,"t"));
    while(ftin >> t0, ftin.good()) tcount++;
    ftin.close();
  }

  open_output(ft,dirsep,"t");
  open_output(fevt,dirsep,"evt");

  if(!restart) {
    remove_dir(Vocabulary->FileName(dirsep,"ekvk"));
    remove_dir(Vocabulary->FileName(dirsep,"transfer"));
  }

  mkdir(Vocabulary->FileName(dirsep,"ekvk"),0xFFFF);
  mkdir(Vocabulary->FileName(dirsep,"transfer"),0xFFFF);

  errno=0;

  if(output)
    open_output(fwk,dirsep,"wk");

  if(movie)
    open_output(fw,dirsep,"w");
}

// wrapper for outcurve routines
class cwrap{
public:
  static Real Spectrum(unsigned int i) {return DNSProblem->Spectrum(i);}
  static Real Dissipation(unsigned int i) {return DNSProblem->Dissipation(i);}
  static Real Pi(unsigned int i) {return DNSProblem->Pi(i);}
  static Real Eta(unsigned int i) {return DNSProblem->Eta(i);}
};

void DNS::Output(int it)
{
  Real E,Z,P;

  w.Set(y);
  ComputeInvariants(w,E,Z,P);
  fevt << t << "\t" << E << "\t" << Z << "\t" << P << endl;

  Complex *y=Y[0];
  if(output) out_curve(fw,y,"w",NY[0]);

  if(movie) OutFrame(it);

  if(spectrum) {
    ostringstream buf;
    Set(T,Y[EK]);
    buf << "ekvk" << dirsep << "t" << tcount;
    open_output(fekvk,dirsep,buf.str().c_str(),0);
    out_curve(fekvk,t,"t");
    out_curve(fekvk,cwrap::Spectrum,"Ek",nshells);
    out_curve(fekvk,cwrap::Dissipation,"nuk*Ek",nshells);
    fekvk.close();
    if(!fekvk) msg(ERROR,"Cannot write to file ekvk");

    Set(T,Y[TRANSFER]);
    buf.str("");
    buf << "transfer" << dirsep << "t" << tcount;
    open_output(ftransfer,dirsep,buf.str().c_str(),0);
    out_curve(ftransfer,t,"t");
    out_curve(ftransfer,cwrap::Pi,"Pi",nshells);
    out_curve(ftransfer,cwrap::Eta,"Eta",nshells);
    ftransfer.close();
    if(!ftransfer) msg(ERROR,"Cannot write to file transfer");
  }

  tcount++;
  ft << t << endl;

  if(rezero && it % rezero == 0 && spectrum) {
    vector2 Y=Integrator->YVector();
    vector T=Y[TRANSFER];
    for(unsigned i=0; i < nshells; i++)
      T[i]=0.0;
    vector S=Y[EK];
    for(unsigned i=0; i < nshells; i++)
      S[i]=0.0;
  }
}
