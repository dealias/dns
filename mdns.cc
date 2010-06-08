#include "options.h"
#include "dns.h"
#include "MultiIntegrator.h"

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
DNS *DNSProblem;

InitialConditionBase *InitialCondition;
ForcingBase *Forcing;

// Global variables for MultiIntegrator.h
MultiProblem *GMProblem; 
unsigned Ngrids;
const char *subintegrator; 

class MDNSVocabulary : public DNSVocabulary {
public:
  const char *Name() {return "MDNS";}
  const char *Abbrev() {return "MDNS";}
  MDNSVocabulary();
  
  //  Table<ForcingBase> *ForcingTable;
  //Table<nltypeBase> *nltypeTable;
  //ForcingBase *NewForcing(const char *& key) {
  //  return ForcingTable->Locate(key);
  //}
};

MDNSVocabulary MDNS_Vocabulary;
DNSVocabulary DNS_Vocabulary;

class MDNS : public MultiProblem {
public:
  MDNS();
  ~MDNS();

  void Initialize();
  void InitialConditions(unsigned Ngrids0);

  void Project(unsigned g);
  void Prolong(unsigned g);

  // class Grid : public DNS {
    // allocate only what we need
    // allocator(?);
  //};
  //array1<Grid *> G;

  void Output(int it);

};

MDNS *MDNSProblem;

MDNS::MDNS() {
  GMProblem=this; // For MultiIntegrator
  grid=0;
  //  exit(1);
  //ExponentialIntegrators(Goy_Vocabulary.IntegratorTable,this);
  //ConservativeIntegrators(Goy_Vocabulary.IntegratorTable,this);
  // FIXME
}
MDNS::~MDNS() {
  // FIXME
}
void MDNS::Initialize() {
  // FIXME
}
void MDNS::InitialConditions(unsigned Ngrids0) {
  // FIXME: this never runs
  MultiProblem::InitialConditions(Ngrids0);
  exit(1);
  if(run=="test")
    run="mtest";
}
void MDNS::Project(unsigned g) {
  // FIXME
}
void MDNS::Prolong(unsigned g) {
  // FIXME
}
void MDNS::Output(int it) {
  // FIXME
}



/****** Vocabulary *****/

MDNSVocabulary::MDNSVocabulary()
{
  Vocabulary=this;
  cout << this << endl;
  //  INTEGRATOR(MultiIntegrator);
  //subintegrator="e_rk3";
  //integrator="multi";
  //  METHOD(MDNS);
}
