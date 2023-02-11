#ifndef __Forcing_h__
#define __Forcing_h__ 1

#define FORCING(key) \
{(void) new Entry<key,ForcingBase>(#key,ForcingTable);}

class ForcingBase {
 public:	
  virtual ~ForcingBase() {}
  virtual const char *Name() {return "None";}
  
  virtual void Init() {}
  virtual void Init(uInt fcount) {}
  virtual bool active(Int i, Int j) {return false;} 
  virtual bool Stochastic(double dt=0.0) {return false;}
  virtual double Force(Complex& w, Complex& S, Int i, Int j) {return 0.0;}
  virtual void ForceMask(Complex& w, Complex& S, Int i, Int j) 
  {
    Force(w,S,i,j);
  }
  virtual double ForceStochastic(Complex& w, Int i, Int j) {return 0.0;}
};

extern ForcingBase *Forcing;
Compare_t ForcingCompare;
KeyCompare_t ForcingKeyCompare;

#endif
