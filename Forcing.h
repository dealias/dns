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
  
  // For stochastic forcing.
  virtual void Force(Array::array2<Complex> &w, const Complex& factor=1.0) {}
  
  // For deterministic forcing.
  virtual void Force(Array::array2<Complex> &w, vector &T) {} 
};

extern ForcingBase *Forcing;
Compare_t ForcingCompare;
KeyCompare_t ForcingKeyCompare;

#endif
