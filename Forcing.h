#ifndef __Forcing_h__
#define __Forcing_h__ 1

#define FORCING(key) \
{(void) new Entry<key,ForcingBase>(#key,ForcingTable);}

enum ForcingType {NONE, CONSTANT, STOCHASTIC};

class ForcingBase {
 public:	
  virtual ~ForcingBase() {}
  virtual const char *Name() {return "None";}
  virtual int Type() {return 0;}
  
  virtual void Set() {}
  virtual void SetStochastic(double dt=0.0) {}
  virtual void Force(Complex& w, double& T, double k) {}
  virtual void ForceStochastic(Complex& w, double& T, double k) {}
};

extern ForcingBase *Forcing;
Compare_t ForcingCompare;
KeyCompare_t ForcingKeyCompare;

#endif
