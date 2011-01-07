#ifndef __rvn_h__
#define __rvn_h__ 1

#include <math.h>
#include "DynVector.h"
#include "Array.h"
using namespace Array;


void check_rvn(DynVector<unsigned> & R2, const unsigned r2, 
	       const unsigned first)
{
  bool found=false;

  unsigned last=R2.Size();
  for(unsigned j=first; j < last; ++j) {
    if(r2 == R2[j]) {
      found=true;
      break;
    }
  }


  /*
  if(R2.Size() > 0) {
    unsigned j = R2.Size();
    do {
      --j;
      if(r2 == R2[j]) {
	found=true;
	j=first; // break;
	break;
      }
    } while(j > first);
  }
  */

  //cout << r2 ;
  if(!found) {
    R2.Push(r2);
    //cout << "pushed" << endl;
  } //else
    //cout << "found" << endl;
}

void findrads(DynVector<unsigned> &R2, array1<unsigned> nr, const unsigned m,
	      const unsigned invis=0)
{
  for(unsigned i=1; i < m; ++i) {
    //unsigned start=nr[(unsigned) floor(sqrt((i-1)/2))];
    double nrstart=floor(sqrt((i-1)/2));
    unsigned start=nr[nrstart > 1 ? (unsigned) nrstart -1 : 0];
    start=0; // FIXME: restore and optimize.
    for(unsigned x=i-1; x <= i; ++x) {
      unsigned x2=x*x;
      unsigned ystopnow= i < x ? i : x;
      for(unsigned y= x == 0 ? 1 : 0; y <= ystopnow; ++y) {
	if(x >= invis || y >= invis) {
	  //cout << "("<<x<<","<<y<<")"<<endl;
	  check_rvn(R2,x2+y*y,start);
	}
      }
    }
    nr[i]=R2.Size();
  }
}
#endif
