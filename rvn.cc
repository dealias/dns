//#include <iostream>
//#include <fstream>
#include "rvn.h"
#include "Array.h"
//using namespace std;
int main()
{
  unsigned N=4;

  DynVector<unsigned> tempR2;
  array1<unsigned> tempnr(N);
  findrads(tempR2,tempnr,N,2);

  cout << tempR2 << endl;

  return 0;
}
