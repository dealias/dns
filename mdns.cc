#include <typeinfo>
#include "dnsbase.h"
#include "MultiIntegrator.h"
#include "Conservative.h"
#include "DynVector.h"
#include "ArrayL.h"

using namespace fftwpp;

const double ProblemVersion=1.0;

#ifndef DEPEND
#if !COMPLEX
#error Navier requires COMPLEX=1 in options.h
#endif
#endif

const char *problem="Multispectral Direct Numerical Simulation of Turbulence";

const char *method="MDNS";
const char *integrator="MultiIntegrator";
const char *ic="Constant";
//const char *linearity="Power";
const char *forcing="WhiteNoiseBanded";

// Vocabulary
int circular=0;
Real nlfactor=0.0;
Real k0=1.0;
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
unsigned rezero=0;
unsigned spectrum=1;
unsigned casimir=0;
Real icalpha=1.0;
Real icbeta=1.0;
enum PRTYPE {NOPR,AREA,POINT};
unsigned prtype=POINT;
int dorescale=false; // TODO: this should be set to true once rescale works
bool overlap=false; // account for overlapping corners of spectrum?
// unused vocab variables to match dns
unsigned movie=0;

// global variables

//***** Global variables for MultiIntegrator.h *****//
MultiProblem *MProblem;
unsigned Ngrids=2;
const char *subintegrator="rk4";

//***** derived loop class *****//
class Mloop : public Hloop {
protected:
  unsigned Invisible;
  Real innerradius, innerradius2;
public:
  void setInvisible(unsigned inv) {Invisible=inv;}

  virtual void setinnerradius(unsigned nold) {
    Real a=max((nold-1.0),0.0);
    innerradius=sqrt(a*a/radix);
    innerradius2=innerradius*innerradius;
    bounds();
  }

  virtual void bounds() {
    Hloop::bounds();
    if(Invisible != 0) {    // subgrids:
      if(circular) {
	Astart=(unsigned) ceil(innerradius)-1;
	Dstart=(unsigned) ceil(innerradius/sqrt2)-1;
      } else { // not circular
	if(radix == 2) {
	  Astart=Invisible/2;
	} else { // radix == 1 or radix == 4
	  Astart=Invisible;
	}
	}
    }
    //cout << Astart << "->" << Astop << endl;
    //cout << Dstart << "->" << Dstop << endl;
  }

  // inner bound for main loop affected by being a subgrid
  virtual unsigned jstart(unsigned I) {
    if(Invisible == 0)
      return 1;
    if(circular) {
      return (unsigned) max(floor(sqrt(innerradius2 - I*I))-1,0.0);
    }
    if(radix == 2)
      return Invisible-I;
    // radix-1 or radix-4
    return I >= Invisible ? 1 : Invisible;
  }
};


//***** Base class for functionoids *****//
class FunctRR {
 public:
  virtual ~FunctRR() {};
  virtual void f(Real w2,Real k2, int n,va_list args) = 0;
};
typedef FunctRR* FunctRRPtr;

// TODO: should these be part of MDNS?
unsigned gN(unsigned N, unsigned g) {
  return radix != 4 ? pow(2,(int) g)*(N+1)-1 : N;
  //  return radix == 1 ? pow(2,(int) g)*(N+1)-1 : N;
};
unsigned gm(unsigned m, unsigned g) {
  return radix != 4 ? pow(2,(int) g)*m : m;
  //return radix == 1 ? pow(2,(int) g)*m : m;
};
unsigned getInvisible(unsigned g) {
  if(g == 0) return 0;
  if(radix==1) return (gN(::Nx,g)+1)/2;
  if(radix==2) return (gN(::Nx,g-1)+1)/2;
  return (gN(::Nx,g-1)+1)/4;  // radix 4
};


//***** MDNS derived from MultiProblem *****//
class MDNS : public MultiProblem {
private:
  unsigned mx, my; // size of data arrays

  unsigned glast;
  void ComputeInvariants(const vector2&, Real&, Real&, Real&);
  array2<Complex> wa, wb; // Vorticity field for projection/prolongation
  array1<unsigned> R2;
  unsigned nshells;
  vector spectra;
  ofstream Mfevt;
  oxstream fekvk;
  ofstream ft;
  oxstream fprolog;
  Real k0;
  int tcount;
  Real etanorm;
public:
  MDNS();
  ~MDNS();

  Real gk(unsigned g) {
    return k0*pow(sqrt((Real) radix),(Real) g);
  };
  void settow(array2<Complex> & w, unsigned g) {
    Set(w,mY[g][OMEGA]);
    G[g]->Dimensiontow(w);
  }

  //***** Grid class based on DNSBase *****//
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
    bool lastgrid;
    unsigned shellsbelow, myshells, mym1;
    DNSBase * smallDNSBase; // for subgrid non-linearity calculation
    unsigned lambda, lambda2; // spacing factor
    
    Mloop gloop;
    Hloop subloop; // for removing doubled interactions

    /*
    bool isvisible(unsigned I, unsigned j){
      return parent->isvisible(I,j,myg);
    }
    */

    /*
    // bounds for spectrum loops
    virtual unsigned diagstart() {
      if(Invisible == 0)
	return 1;
      if(circular) {
	return (unsigned) ceil(((Real) Invisible)/sqrt(2.0));
      }
      if(radix == 2)
	return Invisible/2;
      return Invisible; // radix 1 or 4
    }
    virtual unsigned mainjstart(unsigned I) {
      if(Invisible == 0)
	return 1;
      if(circular) {
	Real a=Invisible-1.0/sqrt(radix);
	return (unsigned) ceil(sqrt(a*a-I*I));
      }
      if(radix == 2)
	return Invisible-I;
      if(I >= Invisible)
	return 1;
      return Invisible;
    }
    virtual unsigned xoriginstart() {
      if(Invisible == 0)
	return 1;
      return Invisible;
    }
    virtual unsigned bottomstart() {
      if(Invisible == 0)
	return 1;
      return Invisible;
    }
    */

  public:
    Grid();
    Grid(unsigned, MDNS *, bool);
    ~Grid();
    void AttachTo(MDNS*, const vector2&);
    void SetParams();

    unsigned myg;

    Real minR2;
    unsigned maxR2;

    array2<Real> tildeB;

    void setcount();
    void setcountoverlap(array1<unsigned>::opt &);
    void setInvisible(unsigned a) {Invisible=a;};
    //unsigned extrashells;

    void settolast() {lastgrid=1;};
    unsigned getnshells() {return nshells;};
    unsigned getNx() {return Nx;};
    unsigned getmy() {return my;};
    unsigned getInvisible() {return Invisible;};
    void setshellsbelow(const unsigned i) {shellsbelow=i;};
    unsigned getshellsbelow() {return shellsbelow;};
    void setmyshells(const unsigned i) {myshells=i;};
    unsigned getmyshells() {return myshells;};
    void Dimensiontow(array2<Complex> & w0) {w0.Dimension(Nx,my);}
    //unsigned nshellbelow;
    unsigned getlambda() {return lambda;}
    unsigned getlambda2() {return lambda2;}

    void InitialConditions(unsigned g);
    //Real gk(Real k, unsigned g) {return k*pow(sqrt((double) radix),g);};
  
    void NonLinearSource(const vector &, const vector &, double);
    void NonConservativeSource(const vector & , const vector &, double);
    void LinearSource(const vector &, const vector &, double);
    void Transfer(const vector2 &, const vector2 &);
    void Spectrum(vector&, const vector&);
    void SpectrumRad2(vector&, const vector&);
    void SpectrumRad4(vector&, const vector&);
    void SpectrumRad4BINNED(vector&, const vector&);
    void SpectrumOverlapRad4BINNED(vector&);
    void Stochastic(const vector2&, double, double);
    
    void killmodes() {exit(1); gloop.killmodes(w);} // FIXME: restore
    void killmodes(Array::array2<Complex> & A) {
      
      //gloop.killmodes(A); // FIXME: restore?
      
      if(circular) { // FIXME: turning this on breaks energy conservation
	for(unsigned i=0; i < Nx; ++i) {
	  //unsigned maxR2=(my-1)*(my-1);
	  unsigned I= i > xorigin ? i-xorigin : xorigin-i;
	  unsigned I2=I*I;
	  vector Ai=A[i];
	  for(unsigned j=0; j < my; ++j) {
	    unsigned r2=I2+j*j;
	    if(r2 >=maxR2) {
	      Ai[j]=0.0;
	    }
	  }
	}
      }

    }

    void ComputeInvariants(const Array::array2<Complex> &,Real&, Real&, Real&);
    void coutw() {cout << "mygrid is "<<myg<<endl<< w << endl;}
    
    void loopwF(const FunctRRPtr,int n,...);
    /*
    class Invariants : public FunctRR {
      void f(const Real w2, const Real k02, int n,va_list args) {
	double * Z=va_arg(args,double *);
	double * E=va_arg(args,double *);
	double * P=va_arg(args,double *);
	*Z += w2;
	*E += w2/k2;
	*P += w2*k2;
      }
    };
    */
  };
  array1<Grid *> G;
  

  // FIXME: remove?
  /*
  bool isvisible(unsigned I, unsigned j, unsigned g) {
    unsigned m=(gN(::Nx,g)+1)/2;
    unsigned Invisible=getInvisible(g);
    if(Invisible==0) {
      if(circular)
	return ((I*I + j*j) <= (m-1)*(m-1));
      return true;
    }
    if(circular) {
      unsigned k2=I*I + j*j;
      Real a=Invisible-1.0; // FIXME: should look more like lamba.asy
      return (k2 >= a*a) && (k2 <= (m-1)*(m-1));
    }
    if(radix == 2)
      return I + j >= Invisible;
    return I >= Invisible || j >= Invisible;
    return true;
  };
  */
  
  void check_rvn(DynVector<unsigned> & R2, const unsigned r2, 
		 const unsigned first)  {
    bool found=false;
    
    unsigned last=R2.Size();
    for(unsigned j=first; j < last; ++j) {
      if(r2 == R2[j]) {
	found=true;
	break;
      }
    }
    if(!found)
      R2.Push(r2);
  }
  
  // FIXME: deprecated
  /*
  void findrads(DynVector<unsigned> &R2, array1<unsigned> nr, unsigned g)  {
    unsigned m=(gN(::Nx,g)+1)/2;
    for(unsigned i=1; i < m; ++i) {
      unsigned start=0; // TODO: calculate better
      for(unsigned x=i-1; x <= i; ++x) {
	unsigned x2=x*x;
	unsigned ystopnow=min(i,x);
	for(unsigned y= x == 0 ? 1 : 0; y <= ystopnow; ++y) {
	  if(isvisible(x,y,g)) {
	    //cout << "("<<x<<","<<y<<")"<<endl;
	    check_rvn(R2,x2+y*y,start);
	  }
	}
      }
      nr[i]=R2.Size();
    }
  }
  */
  
  enum Field {OMEGA,TRANSFER,TRANSFERN,EK,Nfields};
  enum SPEC {NOSPECTRUM, BINNED, INTERPOLATED, RAW}; 
  // copy from DNSBase
  unsigned getnfields(unsigned g) {return Nfields;};
  unsigned getnshells(unsigned g) {
    // this function should only be used during initialization.
    switch(spectrum) {
    case NOSPECTRUM:
      return 0;
      break;
    case BINNED:
      if(radix==1) {
	if(g != glast)
	  return 0;
	else
	  return (unsigned) (hypot(gm(mx,g)-1,gm(my,g)-1)+0.5);
      }

      return gm(my,g)-1;
      break;
    case INTERPOLATED:
      msg(ERROR,"Interpolated spectrum not implemented");
      return 0;
      break;
    case RAW:
      {
	DynVector<unsigned> tempR2;
	array1<unsigned> tempnr(gm(my,g));
	Mloop temploop;
	temploop.setparams(gN(::Nx,g));
	temploop.setInvisible(getInvisible(g));
	temploop.setinnerradius(g == 0 ? 0 : gN(::Nx,g-1));
	temploop.Rloop(tempR2);
	return tempR2.Size();
      }
      break;
    default:
      msg(ERROR,"Invalid spectral choice.");
      return 0;
    }
  };

  //array1<unsigned>::opt count;
  //void setcount();

  Table<InitialConditionBase> *InitialConditionTable;
  void InitialConditions();

  // Source functions
  void Source(const vector2&, const vector2&, double);
  void ConservativeSource(const vector2&, const vector2&, double);
  void NonConservativeSource(const vector2&, const vector2&, double);
  void LinearSource(const vector2&, const vector2&, double);
  void NonLinearSource(const vector2&, const vector2&, double);
  void Transfer(const vector2&, const vector2&);
  void ExponentialSource(const vector2&, const vector2&, double);
  void Stochastic(const vector2&, double, double);
  Real getetanorm() {return etanorm;}
  //void Spectrum(vector&, const vector&);
  //void SpectrumRad4(vector&, const vector&);

  void Project(unsigned ga);
  void Projectrad2point(array2<Complex>&,int,int,array2<Complex>&,int,int,int,			bool overwrite=true);
  void Projectrad1(array2<Complex>&,int,int,array2<Complex>&, int);
  void Prolongrad1(array2<Complex>&,int,int,array2<Complex>&, int);
  void Prolongrad2point(array2<Complex>&,int,int,array2<Complex>&,int,int,int);
  void Projectrad4point(array2<Complex>&,int,int,array2<Complex>&,int,int,int);
  void Prolongrad4point(array2<Complex>&,int,int,array2<Complex>&,int,int,int);
  void Prolong(unsigned gb);
  int Rescale();

  void IndexLimits(unsigned int& start, unsigned int& stop,
		   unsigned int& startT, unsigned int& stopT,
		   unsigned int& startM, unsigned int& stopM) {
    //cout << "grid " << grid << endl;
    unsigned fieldoffset=0;
    for(unsigned g=1; g <= grid; ++g) {
      fieldoffset += getnfields(g-1);
    }
    unsigned offset=grid == 0 ? 0 : Stop(fieldoffset-1);
    start=Start(fieldoffset+OMEGA)-offset;
    stop=Stop(fieldoffset+OMEGA)-offset;
    startT=Start(fieldoffset+TRANSFER)-offset;
    stopT=Stop(fieldoffset+TRANSFER)-offset;
    startM=Start(fieldoffset+EK)-offset;
    stopM=Stop(fieldoffset+EK)-offset;
  }
  
  // Output functions
  Real getSpectrum(unsigned i) {
    switch(spectrum) {
    case BINNED:
      if(radix == 4) {
	for(unsigned g=0; g < Ngrids; ++g) {
	  const unsigned offset=g==0 ? 0 : G[g]->getInvisible()-1;
	  const unsigned firstshell=G[g]->getshellsbelow();
	  const unsigned lastshell=firstshell+G[g]->getmyshells()-1;
	  if (i >=  firstshell &&  i <= lastshell) {
	    unsigned index=i+offset-firstshell;
	    double c=G[g]->count[index];
	    if(c > 0.0)  {
	      return G[g]->T[index].re*twopi/c;
	    }
	    return 0.0; // this should never happen?
	  }
	}
	return 0.0;
      }
      if(radix == 1) {
	// we can just get everything from the most "decimated" grid
	return G[Ngrids-1]->getSpectrum(i);
      }
      msg(ERROR,"Only radix-1 and radix-4 spectra are working right now."); 
      exit(1);
      
      break;
    case INTERPOLATED:
      msg(ERROR,"Interpolated spectrum not done yet.");
      break;
    case RAW:
      {
	// TODO: optimize w. lookup table?
	unsigned k2=R2[i];
	Real val=0.0;
	unsigned c=0;
	for(unsigned g=0; g < Ngrids; ++g) {
	  unsigned gm=G[g]->getmy();
	  unsigned gm2=gm*gm*G[g]->getlambda2();
	  unsigned gnshells=G[g]->getnshells();
	  for(unsigned j=0; j < gnshells ; ++ j) {
	    if(G[g]->getR2(j) == k2) {
	      if(k2 < gm2) {// circular 
		c += G[g]->count[j];
		val += G[g]->T[j].re;
	      }
	    }
	  }
	}
	
	if(c!=0) {
	  return twopi*val/((Real) c);
	}
	return 0.0;
      }
      break;
    default:
      msg(ERROR,"Invalid spectrum choice.");
      break;
    }
    
    msg(ERROR,"Something went really wrong with MDNS::getSpectrum()");
    return 0.0;
  };

  void Computek(DynVector<unsigned>&);
  Real getkb(unsigned i) {
    switch(spectrum) {
    case NOSPECTRUM:
      return 0;
      break;
    case BINNED:
      if(radix==1)
	return k0*i;
      if(radix==4) {
	unsigned lambda=2;
	for(unsigned g=0; g < Ngrids; ++g) {
	  unsigned offset=g==0 ? 0 : G[g]->getInvisible()/lambda;
	  unsigned firstshell=G[g]->getshellsbelow();
	  unsigned lastshell=firstshell+G[g]->getmyshells();
	  if(g==glast) ++lastshell;
	  if (i >=  firstshell &&  i <= lastshell) {
	    unsigned index=i-firstshell+offset;
	    Real kb=gk(g)*(index + 0.5);
	    // 	  cout << "g="<<g<<" k0="<<gk(g) << " index=" << index;
	    // 	  cout << " i=" << i << " kb=" << kb << endl;
	    return kb;
	  }
	}
      }
    case INTERPOLATED:
      msg(ERROR,"Interpolated spectrum not enabled");
      break;
    case RAW:
      if(radix !=4 && radix !=2) 
	msg(ERROR,"only radix-1 and radix-4 cases  with raw spectrum for now.");
      return i == 0 ? 0.5*k0 : k0*sqrt((Real) R2[i-1]);
      break;
    default:
      msg(ERROR,"error in getkb");
      return 0.0; // this should never happen
    }
    msg(ERROR,"error in getkb");
    return 0.0;
  };
  Real getkc(unsigned i) {
    if(spectrum == RAW) 
      return k0*sqrt((Real) R2[i]);
    return 0.5*(getkb(i)+getkb(i+1));
  }

  void Initialize();
  void Output(int it);
  void FinalOutput();
};


//***** Global problem *****//
MDNS *MDNSProblem;


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

class Equipartition : public InitialConditionBase {
public:
  const char *Name() {return "Equipartition";}
  void Set(Complex *w0, unsigned n) {
    unsigned g=MDNSProblem->grid;
    

    unsigned Nx=MDNSProblem->G[g]->getNx();
    unsigned my=MDNSProblem->G[g]->getmy();
    unsigned xorigin=MDNSProblem->G[g]->getxorigin();
    Real k0=MDNSProblem->G[g]->getk0();

    array2<Complex> w(Nx,my,w0);
    w(xorigin,0)=0;
    Real k02=k0*k0;
    for(unsigned i=0; i < Nx; i++) {
      int I=(int) i-(int) xorigin;
      int I2=I*I;
      vector wi=w[i];
      for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
	Real k2=k02*(I2+j*j);
        Real v=icalpha+icbeta*k2;
        v=v ? sqrt(0.5*k2/v) : 0.0;
	wi[j]=Complex(v,v);
      }
    }
    Mloop temploop;
    temploop.setparams(gN(::Nx,g));
    temploop.setInvisible(getInvisible(g));
    temploop.setinnerradius(g == 0 ? 0 : gN(::Nx,g-1));
    temploop.killmodes(w);
    	
  }
};

//***** Forcing *****//
ForcingBase *Forcing;

class None : public ForcingBase {
};

class WhiteNoiseBanded : public ForcingBase {
public:
  const char *Name() {return "White-Noise Banded";}
  void Force(array2<Complex> &w, vector& T, double dt) {
    // FIXME: combine with duplicate code in dns?
    unsigned g=MDNSProblem->grid;
    unsigned Nx=MDNSProblem->G[g]->getNx();
    unsigned my=MDNSProblem->G[g]->getmy();
    unsigned xorigin=MDNSProblem->G[g]->getxorigin();
    Real k02=MDNSProblem->G[g]->getk02();
    Real kmin=max(kforce-0.5*deltaf,0.0);
    Real kmin2=kmin*kmin;
    Real kmax=kforce+0.5*deltaf;
    Real kmax2=kmax*kmax;
    unsigned Invisible=MDNSProblem->G[g]->getInvisible();
    int Invisible2=Invisible*Invisible;

    Real etanorm=MDNSProblem->getetanorm();
    Complex xi=crand_gauss();
    Complex Fk=sqrt(2.0*eta*etanorm);
    Complex fk=Fk*xi;
    double sqrtdt=sqrt(dt);
    Complex diff=sqrtdt*fk;

    
    // TODO: only loop over modes with k in (kmin,kmax)
    //Complex Factor=factor*sqrt(2.0*eta);
    
    for(unsigned i=0; i < Nx; i++) {
      int I=(int) i-(int) xorigin;
      int I2=I*I;
      vector wi=w[i];
      for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
	Real k2=k02*(I2+j*j);
	if(j >= Invisible || I2 >= Invisible2) {
	  if(k2 > kmin2 && k2 < kmax2) {
	    // TODO: enable transfer
	    // T[(unsigned)(sqrt(k2)-0.5)].im += 
	    //   realproduct(Factor,wi[j])+0.5*abs2(Factor);
	    wi[j] += diff;
	  }
	}
      }
    }
  }
};


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
      if(j >= Invisible || I2 >= Invisible2) {
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
  lastgrid=0;
}

MDNS::Grid::Grid(unsigned g, MDNS * parent0, bool islast=0) 
{
  parent=parent0;
  myg=g;
  lastgrid=islast;
  SetParams();
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
  nfields=parent->getnfields(myg);
  w.Set(Y[nfields*myg+OMEGA]);
  Set(T,Y[nfields*myg+EK]);
}

void MDNS::Grid::SetParams()
{
  Nx=gN(::Nx,myg);
  Ny=gN(::Ny,myg);
  k0=parent->gk(myg);
  k02=k0*k0;
  mx=(Nx+1)/2;
  my=(Ny+1)/2;
  mym1=my-1;
  xorigin=mx-1;

  w.Dimension(Nx,my);
  origin=xorigin*my;
  
  Invisible=::getInvisible(myg);
  Invisible2=Invisible*Invisible;
}

void MDNS::Grid::InitialConditions(unsigned g)
{
  myg=g;
  
  // load vocabulary from global variables
  nuH=::nuH;
  nuL=::nuL;
  InitialCondition=MDNS_Vocabulary.NewInitialCondition(ic);

  // check parameters
  if(Nx % 2 == 0 || Ny % 2 == 0) msg(ERROR,"Nx and Ny must be odd");
  if(Nx != Ny) msg(ERROR,"Nx and Ny must be equal");

  //  unsigned Nx0=Nx+xpad;//unused so far...
  //  unsigned Ny0=Ny+ypad; //unused so far...
  //  int my0=Ny0/2+1; //unused so far...

  lambda2=1;
  for(unsigned i=0; i < myg; ++i) 
    lambda2 *= radix;
  lambda=sqrt(lambda2);

  // set up loop structure for this grid
  minR2=0;
  if(myg > 0) {
    unsigned mup=(gN(::Nx,myg-1)+1)/2;
    minR2=(mup-1)*(mup-1)/((Real) radix);
  }

  maxR2=(my-1)*(my-1);
  gloop.setparent(this);
  gloop.setparams(Nx);
  gloop.setInvisible(Invisible);
  gloop.setinnerradius(myg == 0 ? 0 : gN(::Nx,myg-1));
  
  // allocate memory
  cout << "\nGEOMETRY: (" << Nx << " X " << Ny << ")" << endl;
  cout << "\nALLOCATING FFT BUFFERS" << endl;
  block=ComplexAlign(3*Nx*my);
  Convolution=new fftwpp::ImplicitHConvolution2(mx,my,2);
  NLDimension();
  InitialCondition->Set(w,Nx*my);
  fftwpp::HermitianSymmetrizeX(mx,my,xorigin,w);
  gloop.killmodes(w);


  if(myg > 0) { // buffer for calculating small-small-small calculations
    Real sk0=k0;
    if(radix==2) {
      sNx=gN(::Nx,myg-1);
      smx=smy=(sNx+1)/2;
      sk0=k0/sqrt(2.0);
    } else {
      sNx=2*Invisible-1;
      smx=smy=Invisible;
      sxorigin=(sNx+1)/2-1;
    }

    subloop.setparams(sNx);

    unsigned nS=sNx*smy;
      
    xdiff=xorigin-sxorigin;

    wSblock=fftwpp::ComplexAlign(nS);
    SrcwSblock=fftwpp::ComplexAlign(nS);
    wS.Dimension(sNx,smy,wSblock);
    SrcwS.Dimension(sNx,smy,SrcwSblock);
    smallDNSBase=new DNSBase(sNx,smy,sk0); // TODO: share memory block

    if(prtype==AREA) {
      tildeB.Dimension(2*Invisible+1,Invisible+1);
      Allocate(tildeB,(2*Invisible+1)*(Invisible+1));
      for(unsigned i=0; i < 2*Invisible+1; ++i)
	tildeB[i][0]=0.0;
    }
  }

  // set up spectrum parameters
  switch(spectrum) {
  case NOSPECTRUM:
    nshells=0;
    break;
  case BINNED:
    nshells=(unsigned) (hypot(mx-1,my-1)+0.5);
    break;
  case INTERPOLATED:
    msg(ERROR,"Interpolated spectrum not yet enabled.");
    break;
  case RAW:
    {
      DynVector<unsigned> tempR2;
      array1<unsigned> tempnr(my);
      gloop.Rloop(tempR2);
      //parent->findrads(tempR2,tempnr,g);
      //tempR2.sort();
      Allocate(R2,tempR2.Size());
      //Dimension(R2,tempR2.Size());
      for(unsigned i=0; i < R2.Size(); ++i) 
	R2[i]=tempR2[i]*lambda2;

      nshells=R2.Size();
      kval.Dimension(my);
      kval.Allocate(my);
      kval.Load((unsigned) -1);
      for(unsigned I=0; I < mx; ++I) {
	unsigned I2=I*I;
	for(unsigned j=0; j <= I; ++j) {
	  unsigned k2=lambda2*(I2+j*j);
	  for(unsigned a=0; a < nshells; ++a) {
	    if(R2[a] == k2){
	      kval[I][j]=a;
	      //cout << k2 << "," << a << endl;
	    }
	  }
	}
      }
      //cout << "nshells:" << nshells << endl;
      //cout << "R2:\n" << R2 << endl;
      //cout << "kval\n" << kval << endl;
    }
    break;
  default:
    msg(ERROR,"Invalid spectrum type.");
    break;
  }

  // allocate spectral arrays
  if(nshells != 0) {
    nshells=parent->getnshells(myg);
    Dimension(T,nshells);
    Allocate(count,nshells);
    for(unsigned i=0; i < nshells; ++i) {
      T[i]=Complex(0.0,0.0);
      count[i]=0; //FIXME: does this go here?
    }
  }
  setcount();

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

void MDNS::Grid::NonConservativeSource(const vector& Src, const vector& w0, 
				       double)
{
  if(spectrum) {
    vector SrcEK;
    Set(SrcEK,Src);
    Dimension(SrcEK,nshells);
    Spectrum(SrcEK,w0);
  }
}

void MDNS::Grid::Transfer(const vector2 & Src, const vector2 & Y)
{
  // TODO
}

void MDNS::Grid::setcount()
{
  switch(spectrum) {
  case NOSPECTRUM:
    break;
  case BINNED: 
    if(radix == 4) { // FIXME
      //DNSBase::setcountBINNED(); 
    }
    if(lastgrid && radix == 1) { // FIXME
      //DNSBase::setcountBINNED(); 
    }
    if(verbose > 1)  cout << count << endl;
    break;
  case INTERPOLATED: 
    msg(ERROR,"Interpolated spectrum not working right now.");
    break;
  case RAW: { // FIXME
    gloop.Cloop(count,&DNSBase::CountAxes,&DNSBase::CountDiag,
		 &DNSBase::CountMain);
    break;
  }
  default:
    msg(ERROR,"Invalid spectrum choice.");
    break;
  }
}

void pout(int i, int j)
{
  cout << "(" << i << "," << j << ")";
}

void MDNS::Grid::setcountoverlap(array1<unsigned>::opt &Count)
{
  // FIXME: include more spectral choices
  if(spectrum && overlap) {
    Real overlambda=1.0;
    if(radix==4) overlambda=0.5;
    unsigned weight=prtype==AREA ? radix : 1;
    Real kbound=mym1;
    for(unsigned i=0; i < Nx; i++) {
      const int I=(int) i-(int) xorigin;
      const int I2=I*I;
      for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
	if(j >= Invisible || I2 >= Invisible2) {
	  const Real kint=sqrt(I2+j*j);
	  if(kint-0.5 > kbound) {
// 	    cout << "("<<I<<","<<j<<") " << kint << " vs " << kbound
// 		 << " -> " << (unsigned)(overlambda*(kint-0.5)) << endl;
	    Count[(unsigned)(overlambda*(kint-0.5))] += weight;
	  }
	}
      }
    }
  }
}

void MDNS::Grid::Spectrum(vector& SrcEK, const vector& w0)
{
  switch(radix) {
  case 1:
    if(lastgrid) DNSBase::Spectrum(SrcEK,w0);
    break;
  case 2:
    SpectrumRad2(SrcEK,w0);
    break;
  case 4:
    SpectrumRad4(SrcEK,w0);
    break;
  default:
    msg(ERROR,"Invalid radix");
  }
}

void MDNS::Grid::SpectrumRad2(vector& SrcEK, const vector& w0)
{
  switch(spectrum) {
  case NOSPECTRUM:
    break;
  case RAW:
    DNSBase::Spectrum(SrcEK,w0);
    break;
  default:
    msg(ERROR,"Invalid spectrum: only NOSPECTRUM and RAW available");
    exit(1);
    break;
  }
}

void MDNS::Grid::SpectrumRad4(vector& SrcEK, const vector& w0)
{
  switch(spectrum) {
  case NOSPECTRUM:
    break;
  case BINNED:
    SpectrumRad4BINNED(SrcEK,w0);
    if(myg != 0) SpectrumOverlapRad4BINNED(SrcEK);
    break;
  case INTERPOLATED:
    msg(ERROR,"Interpolated spectrum not working yet.");
    break;
  case RAW:
    for(unsigned K=0; K < nshells; K++)    
      SrcEK[K]=Complex(0.0,0.0);
    gloop.Sloop(SrcEK,w,&DNSBase::SpectrumAxes,&DNSBase::SpectrumDiag,
	       &DNSBase::SpectrumMain);
    break;
  default:
    msg(ERROR,"Invalid spectrum");
    exit(1);
    break;
  }
}

void MDNS::Grid::SpectrumRad4BINNED(vector& SrcEK, const vector& w0)
{
  w.Set(w0);
  for(unsigned K=0; K < nshells; K++)
    SrcEK[K]=Complex(0.0,0.0);

  // FIXME: enable different spectral types.
  if(radix != 1) {
    Real kbound=nshells-0.5; // FIXME: does this work?
    //    Real kbound=lastgrid ? hypot(mx,my) : my-0.5;
    for(unsigned i=0; i < Nx; i++) {
      int I=(int) i-(int) xorigin;
      int I2=I*I;
      vector wi=w[i];
      for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
	const Real kint=sqrt(I2+j*j);
	if(kint >= Invisible) {
	  const Real k=k0*kint;
	  const Real w2=abs2(wi[j]);
	  if(kint <= kbound) {
	    SrcEK[(unsigned)(kint-0.5)].re += w2/k;
	    //SrcEK[(unsigned)(k-0.5)].im += nuk(k2)*w2; // TODO
	  }
	}
      }
    }
  } else {
    if(lastgrid) {
      for(unsigned i=0; i < Nx; i++) {
	int I=(int) i-(int) xorigin;
	int I2=I*I;
	vector wi=w[i];
	for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
	  const Real kint=sqrt(I2+j*j);
	  const Real k=k0*kint;
	  const Real w2=abs2(wi[j]);
	  SrcEK[(unsigned)(kint-0.5)].re += w2/k;
	  //SrcEK[(unsigned)(k-0.5)].im += nuk(k2)*w2; // TODO
	}
      }
    }
  }
}

void MDNS::Grid::SpectrumOverlapRad4BINNED(vector& S)
{
  if(spectrum && overlap) {
    MDNSProblem->settow(w,myg);
    Real overlambda=1.0;
    if(radix==4) overlambda=0.5;
    Real weight=prtype==AREA ? radix : 1.0;
    for(unsigned i=0; i < Nx; i++) {
      int I=(int) i-(int) xorigin;
      int I2=I*I;
      vector wi=w[i];
      for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
	if(j >= Invisible || I2 >= Invisible2) {
	  const Real kint=sqrt(I2+j*j);
	  const Real k=k0*kint;
	  const Real overk=1.0/k;
	  const Real w2=abs2(wi[j]);
	  if(kint > my) {
	    S[(unsigned)(overlambda*kint-0.5)].re += weight*w2*overk;
	    //S[(unsigned)(overlambda*kint-0.5)].im += nuk(k2)*w2; 
	  }
	}
      }
    }
  }
}

/****** Vocabulary *****/
MDNSVocabulary::MDNSVocabulary()
{
  Vocabulary=this;

  VOCAB_NOLIMIT(ic,"Initial Condition");
  VOCAB(Nx,1,INT_MAX,"Number of dealiased modes in x direction");
  VOCAB(Ny,1,INT_MAX,"Number of dealiased modes in y direction");
  VOCAB_CONSTANT(movie,0,"Movie flag (off)");
  VOCAB_CONSTANT(casimir,0,"Compute Casimir invariants (off)");
  VOCAB(spectrum,0,3,
	"Spectrum flag (0=off, 1=binned, 2=interpolated, 3=raw)");
  VOCAB(circular,0,1,"kill corner modes?");
  VOCAB(rezero,0,INT_MAX,"Rezero moments every rezero output steps for high accuracy");

  InitialConditionTable=new Table<InitialConditionBase>("initial condition");
  VOCAB(icalpha,0.0,0.0,"initial condition parameter");
  VOCAB(icbeta,0.0,0.0,"initial condition parameter");
  INITIALCONDITION(Zero);
  INITIALCONDITION(Constant);
  INITIALCONDITION(Equipartition); 

  VOCAB(k0,0.0,0.0,"spectral spacing coefficient");
  VOCAB(nuH,0.0,REAL_MAX,"High-wavenumber viscosity");
  VOCAB(nuL,0.0,REAL_MAX,"Low-wavenumber viscosity");
  VOCAB(pH,0,0,"Power of Laplacian for high-wavenumber viscosity");
  VOCAB(pL,0,0,"Power of Laplacian for molecular viscosity");

  VOCAB_NOLIMIT(forcing,"Forcing type");
  ForcingTable=new Table<ForcingBase>("forcing");
  FORCING(None);
  FORCING(WhiteNoiseBanded);
  
  VOCAB(eta,0.0,REAL_MAX,"vorticity injection rate");
  VOCAB(force,(Complex) 0.0, (Complex) 0.0,"constant external force");
  VOCAB(kforce,0.0,REAL_MAX,"forcing wavenumber");
  VOCAB(deltaf,0.0,REAL_MAX,"forcing band width");

  METHOD(MDNS); 

  subintegrator="rk5";
  VOCAB_NOLIMIT(subintegrator,"subintegrator for multi-integrator");
  INTEGRATOR(MultiIntegrator);
  VOCAB(Ngrids,1,INT_MAX,"Number of multispectral grids");
  VOCAB(radix,1,INT_MAX,"Radix number for grid decimation");
  VOCAB(prtype,0,2,"Synchronization scheme (0=none,1=area, 2=point)");
  VOCAB(dorescale,0,1,"Symmetric synchronization? (0=no, 1=hell yes!)");
  VOCAB(nlfactor,0.0,REAL_MAX,"subgrid nonlinear rescaling factor exponent.");
}

MDNS::MDNS() 
{ 
  nfields = Nfields;
  MDNSProblem=this;
  MProblem=this;
  check_compatibility(DEBUG);
  align=sizeof(Complex);
  grid=0;
  ConservativeIntegrators(MDNS_Vocabulary.IntegratorTable,this);
}

MDNS::~MDNS()
{
}

//***** wrappers for output curves *****//
Real curve_Spectrum(unsigned i) {return MDNSProblem->getSpectrum(i);}
Real curve_nuk(unsigned i) {return 0.0;} // TODO
Real curve_kb(unsigned i) {return MDNSProblem->getkb(i);}
Real curve_kc(unsigned i) {return MDNSProblem->getkc(i);}

void MDNS::InitialConditions()
{
  // make sure that the options are cool: 
  if(typeid(*Integrator) != typeid(MultiIntegrator))
    msg(ERROR,"MDNS requires integrator=MultiIntegrator");
  //if(radix==2) msg(ERROR,"radix-2 not implimented");
  if(radix==1 && prtype==AREA) 
    msg(ERROR,"radix-1 only works with pytpe=NONE (0) or POINT (2)");
  
  // Vocabulary
  Ngrids=::Ngrids;
  glast=Ngrids-1;
  saveF=OMEGA;
  MProblem=this;
  k0=::k0;
  Forcing=MDNS_Vocabulary.NewForcing(forcing);
  
  mx=(Nx+1)/2;
  my=(Ny+1)/2;

  // set up MultiProblem
  MultiProblem::InitialConditions(Ngrids);
  for(unsigned g=0; g < Ngrids; ++g) {
    nfields=getnfields(g);
    NY[nfields*g+OMEGA]=gN(Nx,g)*gm(my,g);
    NY[nfields*g+TRANSFER]=0; // getnshells(g); // TODO
    NY[nfields*g+TRANSFERN]=0;
    NY[nfields*g+EK]=getnshells(g);
  }
  //cout << "NY= " << NY << endl;

  // allocate MultiIntegrator and grids
  Allocator(align); 
  G.Allocate(Ngrids);
  for(unsigned g=0; g < Ngrids; ++g) {
    grid=g;
    G[g]=new Grid(g,this,g==glast);
    G[g]->AttachTo(this,Y);
    G[g]->InitialConditions(g);
  }

  // calculate spectrum parameters
  switch(spectrum) {
  case NOSPECTRUM:
    break;
  case BINNED: {
    for(unsigned g=0; g < Ngrids; ++g) {
      G[g]->setcount();
      if(g!=glast) G[g]->setcountoverlap(G[g+1]->count);
    }
    
    Real lambda;
    if(radix==1) lambda=1.0;
    if(radix==4) lambda=2.0;
    unsigned nextra=gm(my,0)-1; // do not include the zero-shell
    nshells = nextra;
    if(verbose >1) cout << "nextra="<<nextra << endl;
    G[0]->setshellsbelow(0);
    G[0]->setmyshells(nextra);
    if(Ngrids > 2) {
      for(unsigned g=1; g < glast; ++g) {
	nextra=gm(my,g) - gm(my,g-1)/lambda;
	G[g]->setshellsbelow(nshells);
	G[g]->setmyshells(nextra);
	nshells += nextra;
	if(verbose >1) cout << "nextra="<<nextra << endl;
      }
    }
    //unsigned gmx=gm(mx,glast); // FIXME: restore?
    //unsigned gmy=gm(my,glast);
    //nextra=(unsigned) (hypot(gmx-1,gmy-1)+0.5*lambda - gm(my,glast-1)/lambda);
    nextra=(unsigned) (getnshells(glast)-G[glast]->getInvisible()+1);
    G[glast]->setshellsbelow(nshells);
    nshells += nextra;
    G[glast]->setmyshells(nextra);
    if(verbose >1) cout << "nextra="<<nextra << endl;
    break;
  }
  case INTERPOLATED:
    msg(ERROR,"Interpolated spectrum not working yet.");
    break;
  case RAW: {
    DynVector<unsigned> tempR2;
    unsigned g=0;
    unsigned gnshells=G[g]->getnshells();
    
    for(unsigned i=0; i < gnshells; ++i) {
      tempR2.Push(G[g]->getR2(i));
    }
        
    for(g=1; g < Ngrids; ++g) {
      gnshells=G[g]->getnshells();
      for(unsigned i=0; i < gnshells; ++i) {
	unsigned temp=G[g]->getR2(i);
	nshells=tempR2.Size();
	bool found=false;
	for(unsigned j=0; j < tempR2.Size(); ++j) {
	  if(temp == tempR2[j]) {
	    found=true;
	    break;
	  }
	}
	if(!found) {
	  tempR2.Push(temp);
	}	    
      }
    }
    nshells=tempR2.Size();
    
    tempR2.sort();
    R2.Allocate(nshells);
    for(unsigned i=0; i < nshells; ++i) 
      R2[i]=tempR2[i];
    break;
  }
  } // switch(spectrum)

  //for (unsigned g=1; g< Ngrids; g++) Project(g);
  // mY not set here for Project to work?


  { // calculate etanorm
    unsigned fcount=0;
    Real kmin=max(kforce-0.5*deltaf,0.0);
    Real kmin2=kmin*kmin;
    Real kmax=kforce+0.5*deltaf;
    Real kmax2=kmax*kmax;

    for(unsigned g=0; g < Ngrids; ++g) {
      unsigned gNx=G[g]->getNx();
      unsigned gmy=G[g]->getmy();
      unsigned gxorigin=G[g]->getxorigin();
      unsigned gInvisible=G[g]->getInvisible();
      Real gk02=G[g]->getk02();

      // FIXME: loop should only be over visible modes.

      for(unsigned i=0; i < gNx; i++) {
	unsigned I= i > gxorigin ? i - gxorigin : gxorigin -1;
	unsigned I2=I*I;
	for(unsigned j=i <= gxorigin ? 1 : 0; j < gmy; ++j) {
	  unsigned k2int=I2+j*j;
	  Real k2=gk02*k2int;
	  if(I >= gInvisible || j >= gInvisible)
	    if(k2 > kmin2 && k2 < kmax2)
	      ++fcount;
	}
      }
      fcount *= 2; // Account for Hermitian conjugate modes.
      if(fcount != 0)
	etanorm=1.0/((Real) fcount);
      else {
	cout << "foring domain not on grid!"<< endl;
	etanorm=1.0;
      }
    }
  }


  // output
  open_output(fprolog,dirsep,"prolog",0);
  out_curve(fprolog,curve_kb,"kb",spectrum ? nshells+1 : 0);
  out_curve(fprolog,curve_kc,"kc",nshells);
  fprolog.close();
  tcount=0;
  //setcount();
  open_output(Mfevt,dirsep,"evt");
  open_output(ft,dirsep,"t");
  if(!restart) remove_dir(Vocabulary->FileName(dirsep,"ekvk"));
  mkdir(Vocabulary->FileName(dirsep,"ekvk"),0xFFFF);
}

void MDNS::Project(unsigned gb)   // FIXME: use Mloop?
{
  if(dorescale==1)
    return;
  if(verbose > 2) cout << "project onto " << G[gb]->myg << endl;
  unsigned ga=gb-1;

  unsigned bInvisible=G[gb]->getInvisible();
  unsigned aNx=G[ga]->getNx();
  unsigned amy=G[ga]->getmy();
  unsigned bNx=G[gb]->getNx();
  unsigned bmy=G[gb]->getmy();
  
  wa.Dimension(aNx,amy);
  wb.Dimension(bNx,bmy);

  Set(wa,mY[ga][OMEGA]);
  Set(wb,mY[gb][OMEGA]);

  G[ga]->killmodes(wa);


  // FIXME: temp
  /*
  wa.Load(1.0);
  wb.Load(2.0);
  G[ga]->killmodes(wa);
  G[gb]->killmodes(wb);
  */

  // debugging code: check if project and prolong conserver energy
#define __checkpr 0
#if __checkpr
  {
    double E=0,Z=0,P=0;
    cout << "\nproject: before: ";
    G[ga]->ComputeInvariants(wa,E,Z,P);
    //cout <<  Z << endl;
    //E=Z=P=0;
    G[gb]->ComputeInvariants(wb,E,Z,P);
    cout << Z << endl;
  }
#endif

  switch(radix) {
  case 1:
    Projectrad1(wa,aNx,amy,wb,bNx);
    break;
  case 2:
    if(prtype == AREA) 
      msg(ERROR,"area-type projection not implemented with radix-2 grids.");
    if(prtype == POINT)
      Projectrad2point(wa,aNx,amy,wb,bNx,bmy,bInvisible);
    break;
  case 4:
    if(prtype==AREA)
      msg(ERROR,"area-type projection no longer implemented");
    if(prtype==POINT)
      Projectrad4point(wa,aNx,amy,wb,bNx,bmy,bInvisible);
    break;
  default:
    msg(ERROR,"Invalid radix.");
    break;
  }

  G[gb]->killmodes(wb);

  //cout << "MDNS::Project\nwa\n:" << wa << "\n wb\n:" << wb ;

#if __checkpr
  {
    double E=0,Z=0,P=0;
    cout << "project: after:  "; 
    G[ga]->ComputeInvariants(wa,E,Z,P);
    //cout << Z << endl;
    //E=Z=P=0;
    G[gb]->ComputeInvariants(wb,E,Z,P);
    cout << Z << endl;
  }
#endif
  //exit(1);
}

void ComplexSet(Complex &a, Complex b) {a=b;}
void ComplexSubfrom(Complex &a, Complex b) {a-=b;}

void MDNS::Projectrad1(array2<Complex>& wa,int aNx,int amy,
		       array2<Complex>& wb, int bNx)
{
  if(prtype != NOPR){
    int axorigin=(aNx-1)/2;
    int bxorigin=(bNx-1)/2;
    int dx=bxorigin-axorigin;
    for(int i=0; i < aNx; i++) {
      vector wai=wa[i];
      vector wbi=wb[i+dx];
      for(int j=0; j < amy; ++j) {
	wbi[j]=wai[j];
      }
    }
  }
}

void MDNS::Projectrad2point(array2<Complex>& wa,int aNx,int amy,
			    array2<Complex>& wb,int bNx,int bmy,int bInvisible,
			    bool overwrite)
{
  int axorigin=(aNx-1)/2;
  int bxorigin=(bNx-1)/2;
  bool iodd=false;
  
  void (*op)(Complex&,Complex)=overwrite ? &ComplexSet : &ComplexSubfrom;

  fftwpp::HermitianSymmetrizeX((aNx+1)/2,amy,axorigin,wa);
  for(int ai=0; ai < aNx; ++ai) {
    int aI=ai-axorigin;
    vector wai;
    Set(wai,wa[ai]);
    Dimension(wai,aNx);
    iodd=!iodd;
    int start=iodd ? 1 : 0;
    if(ai < axorigin && start==0) start=2;
    for(int aj=start; aj < amy; aj+=2) {
      //wai[aj]=Complex(aI,aj);
      int bI=(aI+aj)/2;
      int bj=(aj-aI)/2;
      if(bj >= 0) {
	//cout << "(" << aI+axorigin << "," << aj << ")->"
	//<< "(" << bI+bxorigin << "," << bj << ")" << endl;
	//wb[bI+bxorigin][bj]=Complex(aI,aj);
	//wb[bI+bxorigin][bj]=wai[aj];
	//cout << wai[aj] << " ";
	op(wb[bI+bxorigin][bj],wai[aj]);
	//ComplexSet(wb[bI+bxorigin][bj],wai[aj]);
	//cout << wb[bI+bxorigin][bj] << endl;
      } else {
	//complex conjugate case
	//cout << "(" << aI+axorigin << "," << aj << ")->"
	//<< "(" << bxorigin-bI  << "," << -bj << ")*" << endl;
	//wb[bxorigin-bI][-bj]=Complex(aI,aj); 
	//wb[bxorigin-bI][-bj]=Complex(wai[aj].re,-wai[aj].im);
	//cout <<Complex(wai[aj].re,-wai[aj].im) << " ";
	op(wb[bxorigin-bI][-bj],Complex(wai[aj].re,-wai[aj].im));
	//cout << wb[bxorigin-bI][-bj] << endl;
	//ComplexSet(wb[bxorigin-bI][-bj],wai[aj]);
      }
    }
  }
  fftwpp::HermitianSymmetrizeX((bNx+1)/2,bmy,bxorigin,wb);
  //cout << "wa\n"<<wa << "wb\n"<<wb << endl;
}

void MDNS::Projectrad4point(array2<Complex>& wa,int aNx,int amy,
			    array2<Complex>& wb,int bNx,int bmy,int bInvisible)
{
  int axorigin=(aNx-1)/2;
  int bxorigin=(bNx-1)/2;
  int xstart=bInvisible;
  int xstop=bxorigin+bInvisible;
  int maxR2=(amy-1)*(amy-1); //TODO: get this from the grid?

  int ai=1;
  for(int bi=xstart; bi < xstop; ++bi) {
    vector wai=wa[ai];
    vector wbi=wb[bi];
    int aI=ai-axorigin;
    int aI2=aI*aI;
    for(int j= bi < bxorigin ? 1 : 0 ; j < bInvisible; ++j) {
      if(!circular || (aI2+4*j*j < maxR2))
	wbi[j]=wai[j+j];
    }
    ai += 2;
  }
  fftwpp::HermitianSymmetrizeX((bNx+1)/2,bmy,bxorigin,wb);

  //cout << "MDNS::Projectrad4point\nwa\n:" << wa << "\n wb\n:" << wb ;
}

void MDNS::Prolongrad4point(array2<Complex>& wa,int aNx,int amy,
			    array2<Complex>& wb,int bNx,int bmy,int bInvisible)
{
  int axorigin=(aNx-1)/2;
  int bxorigin=(bNx-1)/2;
  int xstart=bInvisible;
  int xstop=bxorigin+bInvisible;
  int maxR2=(amy-1)*(amy-1);

  int ai=1;
  for(int bi=xstart; bi < xstop; ++bi) {
    int aI=ai-axorigin;
    int aI2=aI*aI;
    vector wai=wa[ai];
    vector wbi=wb[bi];
    for(int j= bi < bxorigin ? 1 : 0 ; j < bInvisible; ++j) {
      if(!circular || (aI2+4*j*j < maxR2))
	wai[j+j]=wbi[j];
    }
    ai += 2;
  }
  fftwpp::HermitianSymmetrizeX((aNx+1)/2,amy,axorigin,wa);

  //cout << "MDNS::Prolongrad4point" << endl;
  //cout << "wa\n:" << wa << "\n wb\n:" << wb 
}

void MDNS::Prolongrad1(array2<Complex>& wa,int aNx,int amy,
		       array2<Complex>& wb,int bNx)
{
  if(prtype != NOPR){
    int axorigin=(aNx-1)/2;
    int bxorigin=(bNx-1)/2;
    int dx=bxorigin-axorigin;
    for(int i=0; i < aNx; i++) {
      vector wai=wa[i];
      vector wbi=wb[i+dx];
      for(int j=0; j < amy; ++j) {
  	wai[j]=wbi[j];
      }
    }
  }
}

void MDNS::Prolong(unsigned ga)
{
  if(dorescale==1)
    return;

  if(verbose > 2) 
    cout << "prolong onto " << G[ga]->myg << endl;
  
  unsigned gb=ga+1;

  unsigned bInvisible=G[gb]->getInvisible();
  unsigned aNx=G[ga]->getNx();
  unsigned amy=G[ga]->getmy();
  unsigned bNx=G[gb]->getNx();
  unsigned bmy=G[gb]->getmy();

  wa.Dimension(aNx,amy);
  wb.Dimension(bNx,bmy);
  Set(wa,mY[ga][OMEGA]);
  Set(wb,mY[gb][OMEGA]);

  G[gb]->killmodes(wb);

  /*
  // FIXME: temp
  wa.Load(1.0);
  wb.Load(2.0);
  G[ga]->killmodes(wa);
  G[gb]->killmodes(wb);
  */

#if __checkpr
  {
    double E=0,Z=0,P=0;
    G[ga]->ComputeInvariants(wa,E,Z,P);
    cout << "prolong: before: ";
    //cout <<  Z << endl;
    //E=Z=P=0;
    G[gb]->ComputeInvariants(wb,E,Z,P);
    cout << Z << endl;
  }
#endif

  switch(radix) {
  case 1:
    Prolongrad1(wa,aNx,amy,wb,bNx);
    break;
  case 2:
    if(prtype == AREA) 
      msg(ERROR,"area-type prolongation not implemented with radix-2 grids.");
    if(prtype == POINT)
      Prolongrad2point(wa,aNx,amy,wb,bNx,bmy,bInvisible);
    break;
  case 4:
    if(prtype==AREA)
      msg(ERROR,"area-type projection no longer implemented");
    if(prtype==POINT)
      Prolongrad4point(wa,aNx,amy,wb,bNx,bmy,bInvisible);
    break;
  default:
    msg(ERROR,"Invalid radix.");
    break;
  }
  G[ga]->killmodes(wa);

#if __checkpr
  {
    double E=0,Z=0,P=0;
    G[ga]->ComputeInvariants(wa,E,Z,P);
    cout << "prolong: after:  ";
    //cout <<  Z << endl;
    //E=Z=P=0;
    G[gb]->ComputeInvariants(wb,E,Z,P);
    cout << Z << endl;
  }  
#endif

  //cout << "MDNS::Prolong\nwa\n:" << wa << "\n wb\n:" << wb ;
}

void MDNS::Prolongrad2point(array2<Complex>& wa,int aNx,int amy,
			    array2<Complex>& wb,int bNx,int bmy,int bInvisible)
{
  //wa.Load(1.0);
  //wb.Load(0.0);
  int axorigin=(aNx-1)/2;
  int bxorigin=(bNx-1)/2;

  fftwpp::HermitianSymmetrizeX((bNx+1)/2,bmy,axorigin,wb);
  int bistop=(int) bNx;
  for(int bi=0; bi < bistop; ++bi) {
    vector wbi;
    Set(wbi,wb[bi]);
    Dimension(wbi,bNx);
    int bI=bi-(int) bxorigin;
    int bjstop=(int)bInvisible-abs(bI);
    for(int bj=bi < (int) bxorigin ? 1 : 0; bj < bjstop; ++bj) {
      int aI=bI-bj;
      int aj=bI+bj;
      if(aj > 0) {
	//cout << "(" << bI+(int)bxorigin << "," << bj  << ")->"
	//<< "(" << aI+(int)axorigin << "," << aj  << ")" << endl;
	//wa[aI+axorigin][aj]=Complex(bI,bj);
	wa[aI+axorigin][aj]=wbi[bj];
      } else {
	//cout << "(" << bI+(int)axorigin << "," << bj  << ")->"
	//<< "(" << (int)axorigin-aI << "," << -aj  << ")*" << endl;
	//complex conjugate case
	//wa[axorigin-aI][-aj]=Complex(bI,bj);
	wa[axorigin-aI][-aj]=Complex(wbi[bj].re,-wbi[bj].im);
	
      }
    }
  }
  fftwpp::HermitianSymmetrizeX((aNx+1)/2,amy,axorigin,wa);
  //cout << "wa\n"<<wa << "wb\n"<<wb << endl;  exit(1);
}

int MDNS::Rescale()
{
  if(dorescale) {
    // FIXME: not done yet
    if(radix == 4 && prtype == POINT) {

      // test code: access current vorticity field
      for(unsigned g=0; g < Ngrids; ++g) {
	array2<Complex> wg;
	settow(wg,g); // FIXME: does settow work?
	wg.Dimension(G[g]->getNx(),G[g]->getmy());
	cout << "wg:" << endl;
	cout << wg << endl;
      }

      // test code: access vorticity field at start of timestep
      for(unsigned g=0; g < Ngrids; ++g) {
	//cout << Ysave[g] << endl; // FIXME: can't access Ysave.
      }

      // test code: find the index for the mode on each grid
      int xstop=8;
      cout << endl;
      //for(int I=0; I <= xstop; I += 2) { // FIXME: restore
      for(int I=0; I <= xstop; I += 1) {
	cout << I << " ";
	int gp=1;
	for(unsigned g=1; g < Ngrids; ++g) {
	  gp *= 2;
	  if((I%gp)==0) {
	    cout << I/gp << " " ;
	  } else { // none of the latter grids are coincident
	    g=Ngrids;
	  }
	}
	cout << endl;
      }
      exit(1);
    
      // find which grids it belongs to.

      // scale it

//       unsigned I=1;

//       for(unsigned i=xstart; i < xstop; ++i) {
// 	//cout << "i="<<i<<" I="<<I<<endl;
// 	vector wai;
// 	Set(wai,wa[I]);
// 	Dimension(wai,aNx);
// 	vector wbi;
// 	Set(wbi,wb[i]);
// 	Dimension(wbi,wb[i]);
// 	for(unsigned j= i < axorigin ? 1 : 0 ; j < bInvisible; ++j) {
// 	  //cout << "j="<<j<<endl;
// 	  wai[j+j]=wbi[j];
// 	}
// 	I += 2;
//       }
//       //cout << "wa:\n"<< wa<< endl << "wb:\n"<< wb << endl;

      return 0;
    }
    return 1; // failure; case not implimented
  }
  return 0; // success
}

void MDNS::Initialize()
{
  Mfevt << "#   t\t\t E\t\t\t Z\t\t\t P" << endl;
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
    out_curve(fekvk,curve_nuk,"nuk*Ek",nshells);
    fekvk.close();
    if(!fekvk) msg(ERROR,"Cannot write to file ekvk");
  }
  tcount++;
  ft << t << endl;
}

void MDNS::ComputeInvariants(const vector2& Y, Real& E, Real& Z, Real& P)
{
  E=Z=P=0.0;
  //cout << "MDNS::ComputeInvariants" << endl;
  for(unsigned g=0; g < Ngrids; ++g) {
    // add up the individual invariants
    Real tempE=0.0, tempZ=0.0, tempP=0.0;

    //G[g]->ComputeInvariants(Y,tempE,tempZ,tempP);

    array2<Complex> w0;
    Set(w0,mY[g][OMEGA]);
    G[g]->Dimensiontow(w0);
    G[g]->ComputeInvariants(w0,tempE,tempZ,tempP);
    Real scale=prtype==AREA ? pow((double) radix,(double) g) : 1.0;
    E += scale*tempE;
    Z += scale*tempZ;
    P += scale*tempP;
  }
}


void MDNS::Grid::ComputeInvariants(const  Array::array2<Complex>& w, 
				   Real& E, Real& Z, Real& P)
{
#if 0
  // FIXME:
  // once project and prolong use Mloop, use Mloop here as well
  gloop.Invariantsloop(w,E,Z,P,&DNSBase::AddInvariants);
#else
  if(circular) {
    //cout << "minR2=" << minR2 << endl;    cout << "maxR2=" << maxR2 << endl;
    //unsigned maxR2=(my-1)*(my-1);
    for(unsigned i=0; i < Nx; i++) {
      int I=(int) i-(int) xorigin;
      int I2=I*I;
      vector wi=w[i];
      for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
	unsigned r2=I2+j*j;
	if((r2 < maxR2) && (r2 >= minR2)) {
	  Real w2=abs2(wi[j]);
	  Real k2=k02*(I2+j*j);
	  Z += w2;
	  E += w2/k2;
	  P += w2*k2;
	}
      }
    }
  } else {
    if(radix == 2) {
      for(unsigned i=0; i < Nx; i++) {
	int I=(int) i-(int) xorigin;
	int I2=I*I;
	vector wi=w[i];
	for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
	  if(abs(I) + j >= Invisible) {
	    Real w2=abs2(wi[j]);
	  Real k2=k02*(I2+j*j);
	  Z += w2;
	  E += w2/k2;
	  P += w2*k2;
	  }
	}
      }
    } else {
      for(unsigned i=0; i < Nx; i++) {
	int I=(int) i-(int) xorigin;
	int I2=I*I;
	vector wi=w[i];
	//cout << "\ni="<<i << ", j=";
	for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
	  if(j >= Invisible || I2 >= Invisible2) {
	    //cout << j<< " ";
	    Real w2=abs2(wi[j]);
	    Real k2=k02*(I2+j*j);
	    Z += w2;
	    E += w2/k2;
	    P += w2*k2;
	  }
	}
      }
    }

  }
  //cout << "grid " << myg << " enstrophy is " << Z << endl;
  //  FunctRRPtr F = new Invariants;
  //  loopwF(F,3,&Z,&E,&P);
#endif
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
  LinearSource(Src,Y,t);
}

/*
void MDNS::setcount() 
{
  count.Allocate(nshells);
  for(unsigned g=0; g < Ngrids; ++g) {
    unsigned gInvisible=G[g]->getInvisible();
    int gInvisible2=gInvisible*gInvisible;
    unsigned gNx=G[g]->getNx();
    unsigned gxorigin=G[g]->getxorigin();
    for(unsigned i=0; i < gNx; i++) {
      int I=(int) i-(int)gxorigin;
      int I2=I*I;
      for(unsigned j=i <= gxorigin ? 1 : 0; j < my; ++j) {
	if(j > gInvisible || I2 > gInvisible2) {
	}
      }
    }
  }
  exit(1);
}
*/

void MDNS::NonConservativeSource(const vector2& Src, const vector2& Y, double t)
{
  vector w0,SrcEK;
  Set(w0,Y[OMEGA]);
  Set(SrcEK,Src[EK]);
  G[grid]->NonConservativeSource(SrcEK,w0,t);
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
  Set(w0,Y[OMEGA]);
  Set(source,Src[OMEGA]);
  G[grid]->LinearSource(source,w0,t);
}

void MDNS::Grid::LinearSource(const vector& source, const vector& w0, double t)
{
  w.Set(w0);
  f0.Set(source);
  for(unsigned i=0; i < Nx; i++) {
    int I=(int) i-(int) xorigin;
    int I2=I*I;
    vector f0i=f0[i];
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      if(j >= Invisible || I2 >= Invisible2) {
	f0i[j] -= nuk(k02*(I2+j*j))*wi[j];
      }
    }
  }
}

void MDNS::NonLinearSource(const vector2& gSrc, const vector2& gY, double t) 
{
  vector w0,source;
  Set(w0,gY[OMEGA]);
  Set(source,gSrc[OMEGA]);
  G[grid]->NonLinearSource(source,w0,t);
}

void MDNS::Grid::NonLinearSource(const vector& wSrc, const vector& wY, double t)
{
  w.Set(wY);
  f0.Set(wSrc);

  //gloop.killmodes(w); // FIXME: restore

  if(circular) { // FIXME: replace with killmodes
    for(unsigned i=0; i < Nx; ++i) {
      //unsigned maxR2=(my-1)*(my-1);
      unsigned I= i > xorigin ? i-xorigin : xorigin-i;
      unsigned I2=I*I;
      vector wi=w[i];
      for(unsigned j=0; j < my; ++j) {
	unsigned r2=I2+j*j;
	if(r2 >=maxR2) {
	  wi[j]=0.0;
	}
      }
    }
  }

  //  cout << "w from Grid::NonLinearSource, myg=" << myg << endl;
  //  cout << w << endl;

  DNSBase::NonLinearSource(f0,w,t);
 
  if(myg > 0) {
    if(radix==2) {
      wS.Load(0.0);
      parent->Prolongrad2point(wS,sNx,smy,w,Nx,my,Invisible);
      smallDNSBase->NonLinearSource(SrcwS,wS,t);
      parent->Projectrad2point(SrcwS,sNx,smy,f0,Nx,my,Invisible,false);
  
    } else { //radix is 1 or 4
      for(unsigned i=0; i < sNx; ++i) {
	vector wi=w[i+xdiff];
	vector wSi=wS[i];
	for(unsigned j=0; j < smy; ++j)
	  wSi[j]=wi[j];
      }
      fftwpp::HermitianSymmetrizeX(smx,smy,sxorigin,wS);
      subloop.killmodes(wS);

      smallDNSBase->NonLinearSource(SrcwS,wS,t);
      subloop.killmodes(SrcwS);
      
      for(unsigned i=0; i < sNx; ++i) {
	vector f0i=f0[i+xdiff];
	vector SrcwSi=SrcwS[i];
	//cout << "i="<< i <<" ,i+xdiff=" <<i+xdiff ;
	for(unsigned j=0; j < smy; ++j)
	  f0i[j] -= SrcwSi[j];
      }
      
    }
    
    if(nlfactor != 0.0) {
      Real fact=pow((double)radix,nlfactor*myg); 
      f0 *= fact;
    }
  } 

  //gloop.killmodes(f0); // FIXME: restore

  if(circular) { // FIXME: replace with killmodes
    for(unsigned i=0; i < Nx; ++i) {
      unsigned I= i > xorigin ? i-xorigin : xorigin-i;
      unsigned I2=I*I;
      vector f0i=f0[i];
      for(unsigned j=0; j < my; ++j) {
	unsigned r2=I2+j*j;
	if(r2 >=maxR2) {
	  f0i[j]=0.0;
	}
      }
    }
  }
  
#if 0
  Real sum=0.0;
  for(unsigned i=0; i < Nx; ++i) {
    unsigned I=i > xorigin ? i-xorigin : xorigin-i;
    unsigned I2=I*I;
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      unsigned r2=I2+j*j;
      if(!(circular) || (r2 < maxR2)) {
	Complex wij=wi[j];
	sum += (f0[i][j]*conj(wij)).re;
      }
    }
  }
  cout << "grid Z conserved by NL? " << sum << endl << endl;
#endif
}

void MDNS::Stochastic(const vector2& gY, double t, double dt) 
{
  G[grid]->Stochastic(gY,t,dt);
}

void MDNS::Grid::Stochastic(const vector2& gY, double t, double dt) 
{
  w.Set(gY[OMEGA]);
  w.Dimension(Nx,my);
  Set(T,gY[TRANSFER]);
  Forcing->Force(w,T,dt);
}

/*
void MDNS::Moments(const vector2& Src, const vector2& Y, double t) 
{
  G[grid]->Moments(Src,Y,t);
}
*/
  
void MDNS::Transfer(const vector2& Src, const vector2& Y) 
{
//  G[grid]->Transfer(Src,Y);
}

