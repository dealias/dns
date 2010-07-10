#include "options.h"
#include "dns.h"
#include "MultiIntegrator.h"

using namespace fftwpp;

const double ProblemVersion=1.0;

#ifndef DEPEND
#if !COMPLEX
#error Navier requires COMPLEX=1 in options.h
#endif
#endif

const char *problem="Multispectral Direct Numerical Simulation of Turbulence";

const char *method="MDNS";
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
unsigned radix=4;
Real eta=0.0;
Complex force=0.0;
Real kforce=1.0;
Real deltaf=1.0;
unsigned movie=0;
unsigned rezero=0;
unsigned spectrum=1;
Real icalpha=1.0;
Real icbeta=1.0;

// global variables for dns.h
int xpad=1;
int ypad=1;

//***** Initial Conditions *****//
InitialConditionBase *InitialCondition;
class Constant : public InitialConditionBase {
public:
  const char *Name() {return "Constant";}
  void Set(Complex *w, unsigned n) {
    for(unsigned i=0; i < n; i++) {
      w[i]=Complex(icalpha,icbeta);
    }
  }
};

//***** Forcing *****//
ForcingBase *Forcing;

//***** Global variables for MultiIntegrator.h *****//
MultiProblem *MProblem;
unsigned Ngrids=2;
const char *subintegrator; 

// ***** Vocabulary ***** //
class MDNSVocabulary : public VocabularyBase {
public:
  const char *Name() {return "MDNS";}
  const char *Abbrev() {return "MDNS";}
  MDNSVocabulary();

  Table<InitialConditionBase> *InitialConditionTable;
  InitialConditionBase *NewInitialCondition(const char *& key) {
    return InitialConditionTable->Locate(key);
  }

  Table<ForcingBase> *ForcingTable;
  ForcingBase *NewForcing(const char *& key) {
    return ForcingTable->Locate(key);
  }
};

MDNSVocabulary MDNS_Vocabulary;


//***** Base class for functionoids *****//
class FunctRR {
 public:
  virtual void f(Real w2,Real k2, int n,va_list args) = 0;
};
typedef FunctRR* FunctRRPtr;

unsigned gN(unsigned N, unsigned g) {
  return radix == 1 ? pow(2,(int) g)*(N+1)-1 : N;
};

unsigned gm(unsigned m, unsigned g) {
  return radix == 1 ? pow(2,(int) g)*m : m;
};


//***** MDNS derived from MultiProblem *****//
class MDNS : public MultiProblem {
private:
  unsigned mx, my; // size of data arrays

  unsigned glast;
  void ComputeInvariants(const vector2 & Y, Real& E, Real& Z, Real& P);
  array2<Complex> wa, wb; // Vorticity field for projection/prolongation
  unsigned nshells;
  vector spectra;
  ofstream Mfevt;
  oxstream fekvk;
  ofstream ft;
  oxstream fprolog;
  int tcount;
  DynVector<unsigned> kvals;
public:

//***** Grid class based on DNS *****//
  class Grid : public DNSBase {
  private:
    unsigned Invisible;
    int Invisible2;
    Complex *wSblock;
    Complex *SrcwSblock;
    array2<Complex> wS, SrcwS;
    unsigned sNx, smx, smy, sxorigin, xdiff;
    void NLDimension();
    MDNS * parent;
    unsigned nfields;

    // for subgrid non-linearity calculation:
    DNSBase * smallDNSBase;

  public:
    Grid();
    Grid(unsigned);
    Grid(MDNS *prob, unsigned g, const vector2& Y0);
    ~Grid();
    void AttachTo(MDNS *prob, const vector2 & Y);
    void SetParams();
    unsigned myg;
  
    unsigned getnshells() {return nshells;};
    unsigned getNx() {return Nx;};
    unsigned getmy() {return my;};
    unsigned getInvisible() {return Invisible;};
    
    void InitialConditions(unsigned g);
    Real gk(Real k, unsigned g) {return k*pow(sqrt((double) radix),g);};
    
    void NonLinearSource(const vector & Src, const vector & Y, double t);
    void LinearSource(const vector & Src, const vector & Y, double t);
    void Transfer(const vector2 & Src, const vector2 & Y);
    void Spectrum(vector& S, const vector& w0);
    
    void ComputeInvariants(const vector2 & Y, Real& E, Real& Z, Real& P);
    void Computek(DynVector<unsigned> & kvals);
    
    void coutw() {cout << "mygrid is"<<myg<<endl<< w << endl;}
  
    void loopwF(const FunctRRPtr,int n,...);
    class Invariants : public FunctRR {
      void f(const Real w2, const Real k02, int n,va_list args) {
	double * Z=va_arg(args,double *);
	double * E=va_arg(args,double *);
	double * P=va_arg(args,double *);
	*Z += w2;
	*E += w2/k02;
	*P += w2*k02;
      }
    };
  };
  array1<Grid *> G;
  
  MDNS();
  ~MDNS();

  enum Field {OMEGA,TRANSFER,EK,Nfields};
  unsigned getnfields(unsigned g) {
    return Nfields;
    //return (g == glast && spectrum) ? Nfields : 1;
  };
  unsigned getNfields() {return Nfields;};
  //unsigned getnshells(unsigned g) {return g == glast ? nshells : 0;};
  unsigned getnshells(unsigned g) {return nshells;}; // FIXME: kludge

  Table<InitialConditionBase> *InitialConditionTable;
  void InitialConditions();

  // Source functions
  void Source(const vector2&, const vector2&, double);
  void ConservativeSource(const vector2& Src, const vector2& Y, double t);
  void NonConservativeSource(const vector2& Src, const vector2& Y, double t);
  void LinearSource(const vector2& Src, const vector2& Y, double t);
  void NonLinearSource(const vector2& Src, const vector2& Y, double t);
  void Transfer(const vector2& Src, const vector2& Y);
  void ExponentialSource(const vector2& Src, const vector2& Y, double t);
  void Stochastic(const vector2& Y, double t, double dt);
  
  void Project(unsigned ga);
  void Prolong(unsigned gb);
  
  // Output functions
  Real getSpectrum(unsigned i) {return Y[EK][i].re;}; 
  void Computek(DynVector<unsigned>&);
  Real getk(unsigned i) {return (Real) i;};

  void FinalOutput();
  void Initialize();
  void Output(int it);
};

//***** Global problem *****//
MDNS *MDNSProblem;

// ***** Grid functions *****//
void MDNS::Grid::loopwF(const FunctRRPtr F,int n,...)
{
  va_list args;
  va_start(args,n);
  for(unsigned i=0; i < Nx; i++) {
    int I=(int) i-(int) xorigin;
    int I2=I*I;
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      if(j > Invisible || I2 > Invisible2) {
	Real w2=abs2(wi[j]);
	Real k2=k02*(I2+j*j);
	F->f(w2,k2,n,args);
      }
    }
  }
  va_end(args);
};


MDNS::Grid::Grid()
{
}

MDNS::Grid::Grid(unsigned g) 
{
  myg=g;
  SetParams();
}

MDNS::Grid::Grid(MDNS *prob, unsigned g, const vector2 & Y)
{
  parent=prob;
  //  w.Set(Y[OMEGA]);
  w.Set(Y[nfields*myg+OMEGA]); // FIXME: is this right?
}

MDNS::Grid::~Grid()
{
  fftwpp::deleteAlign(block);
  if(myg > 0) {
    fftwpp::deleteAlign(wSblock);
    fftwpp::deleteAlign(SrcwSblock);
  }
}

void MDNS::Grid::AttachTo(MDNS *prob, const vector2 & Y) 
{
  parent=prob;
  nfields=parent->getNfields();
  w.Set(Y[nfields*myg+OMEGA]); // FIXME: is this right?
}

void MDNS::Grid::SetParams()
{
  Nx=gN(::Nx,myg);
  Ny=gN(::Ny,myg);
  k0=gk(1.0,myg); // assumes k0 is set to one for undecimated grid
  k02=k0*k0;
  mx=(Nx+1)/2;
  my=(Ny+1)/2;
  xorigin=mx-1;

  w.Dimension(Nx,my);
  origin=xorigin*my;
  
  if(myg == 0) {
    Invisible=Invisible2=0;
  } else {
    Invisible=(gm(::Nx,myg-1)/2)*gk(1.0,myg-1)/gk(1.0,myg);
    Invisible2=Invisible*Invisible;
  }
}

void MDNS::Grid::InitialConditions(unsigned g)
{
  // load vocabulary from global variables
  myg=g;
  xpad=::xpad;
  ypad=::ypad;
  nuH=::nuH;
  nuL=::nuL;
  pH=::pH;
  pL=::pL;
  spectrum=::spectrum;
  //  SetParams();

  if(Nx % 2 == 0 || Ny % 2 == 0) msg(ERROR,"Nx and Ny must be odd");
  if(Nx != Ny) msg(ERROR,"Nx and Ny must be equal");
  if(radix != 1 && radix != 4) 
    msg(ERROR,"only radix 1 (trivial) or 4 decimations enabled");

  nshells=MDNSProblem->getnshells(myg);

  cout << "\nGEOMETRY: (" << Nx << " X " << Ny << ")" << endl;

  cout << "\nALLOCATING FFT BUFFERS" << endl;

  Dimension(T,nshells);

  //  unsigned Nx0=Nx+xpad;//unused so far...
  //  unsigned Ny0=Ny+ypad; //unused so far...
  //  int my0=Ny0/2+1; //unused so far...


  block=ComplexAlign(3*Nx*my);
  Convolution=new fftwpp::ImplicitHConvolution2(mx,my,2);
  NLDimension();

  if(myg > 0) { 
    // buffer for calculating small-small-small calculations
    sNx=Invisible;
    smx=smy=(Invisible+1)/2;
    sxorigin=(sNx+1)/2-1;
    unsigned nS=sNx*smy;
    xdiff=xorigin-sxorigin;

    wSblock=fftwpp::ComplexAlign(nS);
    SrcwSblock=fftwpp::ComplexAlign(nS);
    wS.Dimension(sNx,smy,wSblock);
    SrcwS.Dimension(sNx,smy,SrcwSblock);
    smallDNSBase=new DNSBase(sNx,smy,k0); // TODO: share memory block
  }

  InitialCondition=MDNS_Vocabulary.NewInitialCondition(ic);
  InitialCondition->Set(w,Nx*my);

  fftwpp::HermitianSymmetrizeX(mx,my,xorigin,w);
}

void MDNS::Grid::NLDimension()
{
  w.Dimension(Nx,my);
  f0.Dimension(Nx,my);
  f1.Dimension(Nx,my,block);
  g0.Dimension(Nx,my,block+Nx*my);
  g1.Dimension(Nx,my,block+2*Nx*my);
  F[0]=f0;
  F[1]=f1;
  G[0]=g0;
  G[1]=g1;
}

void MDNS::Grid::LinearSource(const vector& source, const vector& w0, double t)
{
  // FIXME: overcounts
  DNSBase::LinearSource(source,w0,t);
}

void MDNS::Grid::NonLinearSource(const vector& source, const vector& w0,
				 double t)
{
  w.Set(w0);
  f0.Set(source);

  DNSBase::NonLinearSource(source,w0,t);

  if(myg > 0) {
    //copy overlapping modes to wS
    for(unsigned i=0; i < sNx; ++i) {
      vector wi=w[i+xdiff];
      vector wSi=wS[i];
      for(unsigned j=0; j < smy; ++j)
	wSi[j]=wi[j];
    }
    fftwpp::HermitianSymmetrizeX(smx,smy,sxorigin,wS);

    // find the nonlinear interaction for small-small-small
    smallDNSBase->NonLinearSource(SrcwS,wS,t); 
 
    // subtract small-small-small source
    for(unsigned i=0; i < sNx; ++i) {
      vector f0i=f0[i+xdiff];
      vector SrcwSi=SrcwS[i];
      for(unsigned j=i <= xorigin ? 1 : 0; j < smy; ++j)
	f0i[j] -= SrcwSi[j];
    }
    fftwpp::HermitianSymmetrizeX(mx,my,xorigin,f0);
  }
}

void MDNS::Grid::Transfer(const vector2 & Src, const vector2 & Y)
{
  // FIXME
}

void MDNS::Grid::Spectrum(vector& S, const vector& w0)
{
  w.Set(w0);
  
  // FIXME: loop over only visible modes!

  // Compute instantaneous angular sum over each circular shell.
  for(unsigned i=0; i < Nx; i++) {
    int I=(int) i-(int) xorigin;
    int I2=I*I;
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      Real k2=k02*(I2+j*j);
      Real k=sqrt(k2);
      Real w2=abs2(wi[j]);
      S[(unsigned)(k-0.5)] += Complex(w2/k,nuk(k2)*w2); //nuk set?
    }
  }
}

void MDNS::Grid::Computek(DynVector<unsigned> & kvals)
{
  
  /*
  for(unsigned i=xorigin; i < Nx; ++i) {
    int I2=(int) (i - xorigin)*(i - xorigin);
    for(unsigned j=0; j < my; ++j) { // FIXME: include origin?
      Real k=sqrt(k02*(I2+j*j));
      bool set=false;
      for(unsigned a=0; a < kvals.Size() && !set; ++a) {
	set=false;
	if(k==kvals[a])
	  set=true;
      }
      if(!set)
	kvals.Push(k);
    }
  }
  cout << kvals.Size() << endl;
  */
}

void MDNS::Grid::ComputeInvariants(const vector2 & Y, Real& E, Real& Z, Real& P)
{
  w.Set(Y[nfields*myg+OMEGA]); // FIXME: is this right?
  //  w.Set(Y[OMEGA]);
  //Dimension(w,Y[OMEGA]); // FIXME: breaks makeopt
  FunctRRPtr F = new Invariants;
  loopwF(F,3,&Z,&E,&P);
}

/****** Vocabulary *****/
MDNSVocabulary::MDNSVocabulary()
{
  Vocabulary=this;

  VOCAB_NOLIMIT(ic,"Initial Condition");
  VOCAB(Nx,1,INT_MAX,"Number of dealiased modes in x direction");
  VOCAB(Ny,1,INT_MAX,"Number of dealiased modes in y direction");
  //  VOCAB(movie,0,1,"Movie flag (0=off, 1=on)");
  VOCAB(spectrum,0,1,"Spectrum flag (0=off, 1=on)");
  VOCAB(rezero,0,INT_MAX,"Rezero moments every rezero output steps for high accuracy");

  InitialConditionTable=new Table<InitialConditionBase>("initial condition");
  VOCAB(icalpha,0.0,0.0,"initial condition parameter");
  VOCAB(icbeta,0.0,0.0,"initial condition parameter");
  INITIALCONDITION(Zero);
  INITIALCONDITION(Constant);
  //  INITIALCONDITION(Equipartition); 

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
  //FORCING(None);
  //FORCING(WhiteNoiseBanded);

  METHOD(MDNS); 

  subintegrator="rk5";
  VOCAB_NOLIMIT(subintegrator,"subintegrator for multi-integrator");
  INTEGRATOR(MultiIntegrator);
  VOCAB(Ngrids,1,INT_MAX,"Number of multispectral grids");
  VOCAB(radix,1,INT_MAX,"Radix number for grid decimation");
}

MDNS::MDNS() 
{
  nfields = Nfields;
  MDNSProblem=this;
  MProblem=this;
  check_compatibility(DEBUG);
  align=sizeof(Complex);
}

MDNS::~MDNS() 
{
}

//***** wrappers for output curves *****//
Real curve_Spectrum(unsigned i) {return MDNSProblem->getSpectrum(i);}
Real curve_k(unsigned i) {return MDNSProblem->getk(i);}

void MDNS::InitialConditions()
{
  //***** Vocabulary *****//
  Ngrids=::Ngrids;
  glast=Ngrids-1;
  saveF=OMEGA;
  MultiProblem::InitialConditions(Ngrids);
  MProblem=this;

  mx=(Nx+1)/2;
  my=(Ny+1)/2;

  G.Allocate(Ngrids);
  for(unsigned int g=0; g < Ngrids; ++g) {
    G[g]=new Grid(g);
  }

  if(spectrum) { 
    //Computek(kvals);
    //nshells=kvals.Size();
    nshells=spectrum ? 
      (unsigned) (G[glast]->getk0()*hypot(gm(mx,glast)-1,gm(my,glast)-1)+0.5) 
      : 0; // FIXME: k0??
  } else
    nshells=0;

  Allocate(spectra,nshells);
  for (unsigned i=0; i < nshells; ++i) 
    spectra[i]=0.0;

  for(unsigned g=0; g < Ngrids; ++g) {
    nfields=getNfields();
    NY[nfields*g+OMEGA]=gN(Nx,g)*gm(my,g);
    NY[nfields*g+TRANSFER]=0; // getnshells(g); FIXME 
    NY[nfields*g+EK]=getnshells(g);
  }

  Allocator(align); // allocate MultiIntegrator

  for(unsigned int g=0; g < Ngrids; ++g)
    G[g]->AttachTo(this,Y);

  for(unsigned g=0; g < Ngrids; ++g)
    G[g]->InitialConditions(g);

  // FIXME: kludge
  for(unsigned g=0; g < Ngrids; ++g) {
    unsigned stop=G[g]->getNx()*G[g]->getmy();
    for(unsigned i=0; i < stop; ++i)
      Y[nfields*g+OMEGA][i]=Complex(icalpha,icbeta);
  }

  //***** output *****//
  // FIXME: is anyone else writing to fprolog?
  open_output(fprolog,dirsep,"prolog",0);
  out_curve(fprolog,curve_k,"k",G[glast]->getnshells());
  fprolog.close();

  tcount=0;
  open_output(Mfevt,dirsep,"evt");
  open_output(ft,dirsep,"t");
  if(!restart) remove_dir(Vocabulary->FileName(dirsep,"ekvk"));
  mkdir(Vocabulary->FileName(dirsep,"ekvk"),0xFFFF);
}

void MDNS::Project(unsigned gb) 
{
  return; // FIXME: temp
  //  cout << "project onto " << G[gb]->myg << endl;
  unsigned ga=gb-1;

  unsigned  aInvisible=G[ga]->getInvisible();
  unsigned  bInvisible=G[gb]->getInvisible();
  unsigned  axorigin=G[ga]->getxorigin();
  unsigned  bxorigin=G[gb]->getxorigin();
  unsigned aNx=G[ga]->getNx();
  unsigned amy=G[ga]->getmy();
  unsigned  amx=G[ga]->getmx();
  unsigned  bmx=G[gb]->getmx();
  unsigned bNx=G[gb]->getNx();
  unsigned bmy=G[gb]->getmy();
  unsigned dx=bxorigin-axorigin;
  Real ak02=G[ga]->getk02();
  Real bk02=G[gb]->getk02();

  wa.Dimension(aNx,amy);
  wb.Dimension(bNx,bmy);
  Set(wa,mY[ga][OMEGA]);
  Set(wb,mY[gb][OMEGA]);
  
  const unsigned xstart=bxorigin-bInvisible;
  const unsigned xstop=bxorigin+bInvisible;

  if(radix == 1) {
    for(unsigned int i=xstart; i <= xstop; i++) {
      vector wai=wa[i-dx];
      vector wbi=wb[i];
      for(unsigned int j=i <= bxorigin ? 1 : 0; j <= bInvisible; ++j) {
  	wbi[j]=wai[j];
      }
    }
  }

  if(radix == 4) {
    for(unsigned i=xstart; i <= xstop; ++ i) {
      int I=(int)i - (int)bxorigin;
      vector wai=wa[2*I+axorigin];
      vector waim=wa[2*I+axorigin-1];;
      vector waip=wa[2*I+axorigin+1];
      vector wbi=wb[i];
      
      for(unsigned j= i <= axorigin ? 1 : 0 ; j <= bInvisible; ++j) {
	
 	const int aJ=2*j;
	const int aJp=aJ+1;
	const int aJm=aJ==0? aJp : aJ-1;

	// project co-incident point
	wai[aJ]=Complex(2*I,aJ);
	
	//project points on same row/column
	waim[aJ]=Complex(2*I-1,aJ);
	waip[aJ]=Complex(2*I+1,aJ);
	if(aJ) {// Hermiticity case. 
	  wai[aJm]=Complex(2*I,aJm); 
	} else {
	  wa[axorigin-2*I][aJp]=Complex(-2*I,aJm); //--
	}
	wai[aJp]=Complex(2*I,aJp);
	
	//project points on diagnols.
	waim[aJp]=Complex(2*I-1,aJp); // -+
	waip[aJp]=Complex(2*I+1,aJp); // ++
	if(aJ) {// Hermiticity case.
	  waim[aJm]=Complex(2*I-1,aJm);
	  waip[aJm]=Complex(2*I+1,aJm);
	} else {
	  wa[axorigin-2*I+1][aJm]=Complex(-2*I+1,aJm); // --
	  wa[axorigin-2*I-1][aJm]=Complex(-2*I-1,aJm); // +-
	}
	
      }
    }

    for(unsigned i=0; i < aNx; i += 2) {
      int I= (int) i - (int) bxorigin;
      int bI=I/2;
      unsigned bi=bI + bxorigin;
      for(unsigned j= i < axorigin? 2 : 0; j <  amy; j +=2) {
	unsigned bj=j/2;
	Real bZ=abs2(wb[bi][bj]);
      }
    }
  }

  

  HermitianSymmetrizeX(bNx,bmy,bxorigin,wb); 
}

void MDNS::Prolong(unsigned ga)
{
  return;  // FIXME: temp

  //  cout << "prolong onto " << G[ga]->myg << endl;
  unsigned gb=ga+1;
  
  int aInvisible=(int) G[ga]->getInvisible();
  int axorigin=(int) G[ga]->getxorigin();
  int amx=(int) G[ga]->getmx();
  int dx=(int) G[gb]->getxorigin()-axorigin;

  wa.Dimension(G[ga]->getNx(),G[ga]->getmy());
  wb.Dimension(G[gb]->getNx(),G[gb]->getmy());

  // FIXME: code doesn't work under both make and makeopt.
  // for make:
  ///*
  array1<array1<Complex> > va=mY[ga];
  Set(wa,va[OMEGA]);
  array1<array1<Complex> > vb=mY[gb];
  Set(wb,vb[OMEGA]);
  //*/

  // for makeopt:
  /*
  array1<Complex *> va=mY[ga];
  wa.Set(va[OMEGA]);
  array1<Complex *> vb=mY[gb];
  wb.Set(vb[OMEGA]);
  */



  if(radix == 1) {
    const int xstart=amx-aInvisible;
    const int xstop=amx-aInvisible;
    for(int i=xstart; i < xstop; i++) {
      vector wai=wa[i];
      vector wbi=wb[i+dx];
      for(int j=i <= axorigin ? 1 : 0; j < aInvisible; ++j) {
	wai[j]=wbi[j];
      }
    }
  }
  if(radix == 4) {
    // TODO: radix4
  }

  fftwpp::HermitianSymmetrizeX(amx,G[ga]->getmy(),axorigin,wa);
  // maybe only on the overlapping modes?

  // FIXME: copy stuff from spectra onto just onto lastgrid's Src[EK]
}

void MDNS::Initialize()
{
  Mfevt << "#   t\t\t E\t\t\t Z" << endl;
}

void MDNS::Output(int it)
{
  Real E,Z,P;

  ComputeInvariants(Y,E,Z,P);
  Mfevt << t << "\t" << E << "\t" << Z << "\t" << P << endl;
    
  if(spectrum) {
    ostringstream buf;
    buf << "ekvk" << dirsep << "t" << tcount;
    open_output(fekvk,dirsep,buf.str().c_str(),0);
    out_curve(fekvk,t,"t");
    out_curve(fekvk,curve_Spectrum,"Ek",nshells);
    out_curve(fekvk,curve_Spectrum,"nuk*Ek",nshells); // FIXME: Temp
    fekvk.close();
    if(!fekvk) msg(ERROR,"Cannot write to file ekvk");
  }
  tcount++;
  ft << t << endl;
}

void MDNS::ComputeInvariants(const vector2& Y, Real& E, Real& Z, Real& P)
{
  E=Z=P=0.0;

  Real tempE=0.0, tempZ=0.0, tempP=0.0;
  for(unsigned g=0; g < Ngrids; ++g) {
    // add up the individual invariants
    G[g]->ComputeInvariants(Y,tempE,tempZ,tempP);
    Real scale=pow((double) radix,(double) g);
    E += scale*tempE;
    Z += scale*tempZ;
    P += scale*tempP;
  }
}

void MDNS::Computek(DynVector<unsigned> & kvals) 
{
  /*
  for(unsigned g=0; g < Ngrids; ++g)
    G[g]->Computek(kvals);

  // sort kvals with quick-'n'-dirty bubble sort:
  bool flag=1;
  unsigned temp;
  unsigned n=kvals.Size();
  for(unsigned i=1; (i < n) && flag; i++) {
    flag=0;
    for(unsigned j=0; j < (n-1); j++) {
      if(kvals[j+1] < kvals[j]) { 
	temp= kvals[j];    
	kvals[j] = kvals[j+1];
	kvals[j+1] = temp;
	flag=1;
      }
    }
  }
  */
}

void MDNS::FinalOutput()
{
  Real E,Z,P;
  ComputeInvariants(Y,E,Z,P);
  cout << endl;
  cout << "Energy = " << E << newl;
  cout << "Enstrophy = " << Z << newl;
  cout << "Palenstrophy = " << P << newl;
}

void MDNS::Source(const vector2& Src, const vector2& Y, double t)
{
  ConservativeSource(Src,Y,t);
  NonConservativeSource(Src,Y,t);
}

void MDNS::ConservativeSource(const vector2& Src, const vector2& Y, double t) 
{
  NonLinearSource(Src,Y,t);
  Transfer(Src,Y);
  LinearSource(Src,Y,t); // is this really a conservative source?
}

void MDNS::NonConservativeSource(const vector2& Src, const vector2& Y, double t)
{
  if(spectrum) {
    for(unsigned g=0; g < Ngrids; ++g) {
      if(g==0) { // zero the spectrum
	for(unsigned i=0; i < nshells; ++i)
	  spectra[i]=0;
      }
      vector w0,source;
      Dimension(w0,Y[OMEGA]);
      G[g]->Spectrum(spectra,w0);
      
      if(g==glast) {
	vector S;
	Dimension(S,Src[EK]);
	for (unsigned i=0; i < nshells; ++i)
	  S[i]=spectra[i];  // might this instead be a swap?
      }
    }
  }
  //Moments(Src,Y,t);
}

void MDNS::ExponentialSource(const vector2& Src, const vector2& Y, double t) 
{
  NonLinearSource(Src,Y,t);
  Transfer(Src,Y);
  NonConservativeSource(Src,Y,t);
}

void MDNS::LinearSource(const vector2& Src, const vector2& Y, double t) 
{
  vector w0,source;
  Dimension(w0,Y[OMEGA]);
  Dimension(source,Src[OMEGA]);
  G[grid]->LinearSource(source,w0,t);
}

void MDNS::NonLinearSource(const vector2& Src, const vector2& Y, double t) 
{
  vector w0,source;
  Dimension(w0,Y[OMEGA]);
  Dimension(source,Src[OMEGA]);
  G[grid]->NonLinearSource(source,w0,t);
}

void MDNS::Stochastic(const vector2& Y, double t, double dt) 
{
  //G[grid]->Stochastic(T,y,t,dt);
}

/*
void MDNS::Moments(const vector2& Src, const vector2& Y, double t) 
{
  G[grid]->Moments(Src,Y,t);
}
*/
  
void MDNS::Transfer(const vector2& Src, const vector2& Y) 
{
  G[grid]->Transfer(Src,Y);
}

