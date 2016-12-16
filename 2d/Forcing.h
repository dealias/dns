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
  virtual bool active(int i, int j) {return false;} 
  virtual bool Stochastic(double dt=0.0) {return false;}
  virtual void Force(Complex& w, Complex& S, double& T, int i, int j) {}
  virtual void ForceMask(Complex& w, Complex& S, double& T, int i, int j) 
  {
    Force(w,S,T,i,j);
  }
  virtual void ForceStochastic(Complex& w, double& T, int i, int j) {}
};

extern ForcingBase *Forcing;
Compare_t ForcingCompare;
KeyCompare_t ForcingKeyCompare;

#endif
