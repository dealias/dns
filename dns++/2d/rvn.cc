#include <iostream>
#include <fstream>
#include "DynVector.h"
#include <math.h>

using namespace std;

void check(DynVector<int> & R2, const int r2, const int first)
{
  bool found=false;
  int last=R2.Size();
  for(int j=last > 0 ? last-1 : last; j >= first ; --j) {
    //cout << "j="<<j << endl;
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
  int N=400;

  ofstream file;
  file.open ("rvn.dat");

  DynVector<int> R2;
  DynVector<int> nr(N);
  
  for(int i=1; i < N; ++i) {
    cout << "i="<< i << endl;

    int start=nr[(int) floor(sqrt((i-1)/2))];

    for(int x=i-1; x <= i; ++x) {
      int x2=x*x;
      int ystopnow=min(i,x);
      for(int y=0; y <= ystopnow; ++y) {
	check(R2,x2+y*y,start);
      }
    }

    //cout << R2 << endl;
    file << i << "\t" << R2.Size() << "\n";
    nr[i]=R2.Size();
    cout << "nr=" << R2.Size() << endl << endl;
    
  }
  file.close();
  cout << nr.Size();

  return 0;
}
