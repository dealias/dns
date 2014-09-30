#ifndef __Forcing_h__
#define __Forcing_h__ 1

#define FORCING(key) \
{(void) new Entry<key,ForcingBase>(#key,ForcingTable);}

class ForcingBase {
 public:	
  virtual ~ForcingBase() {}
  virtual const char *Name() {return "None";}
  
  virtual void Init() {}
  virtual void Init(unsigned fcount) {}
  virtual bool active(double k) {return false;} 
  virtual bool Stochastic(double dt=0.0) {return false;}
  virtual void Force(Complex& w, double& T, double k) {}
  virtual void ForceStochastic(Complex& w, double& T, double k) {}
};

extern ForcingBase *Forcing;
Compare_t ForcingCompare;
KeyCompare_t ForcingKeyCompare;

#endif
