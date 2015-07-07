#ifndef __InitialCondition_h__
#define __InitialCondition_h__ 1

#define INITIALCONDITION(key) \
{(void) new Entry<key,InitialConditionBase>(#key,InitialConditionTable);}

class InitialConditionBase {
 public:	
  virtual ~InitialConditionBase() {};
  virtual const char *Name() {return "None";}
  virtual Var Value(Real k) {return 0;}
};

extern InitialConditionBase *InitialCondition;
Compare_t InitialConditionCompare;
KeyCompare_t InitialConditionKeyCompare;

#endif
