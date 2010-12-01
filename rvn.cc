#include <iostream>
#include <fstream>
#include "DynVector.h"
using namespace std;

void check(DynVector<int> & R2, const int r2)
{
  bool found=false;
  int last=R2.Size();
  for(int j=last-1; j >= 0 ; --j) {
    if(R2[j]==r2) {
      found=true;
      break;
    }
  }

  if(!found)
    R2.Push(r2);
}

int main()
{
  unsigned N=1000;

  ofstream file;
  file.open ("rvt.txt");

  DynVector<int> R2;
  
  for(int i=1; i < N; ++i) {
    cout << "i="<< i << endl;
    for(int x=i-1; x <= i; ++x) {
      int x2=x*x;
      for(int y=0; y <= x; ++y) {
	check(R2,x2+y*y);
      }
    }
    for(int x=i-1; x <= i; ++x) {
      int x2=x*x;
      for(int y=i-1; y <= i; ++y) {
	check(R2,x2+y*y);
      }
    }

    file << i << "\t" << R2.Size() << "\n";
    cout << "nr=" << R2.Size() << endl << endl;
    
  }

  file.close();

  return 0;
}

