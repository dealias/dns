#include "Array.h"
#include "DynVector.h"
using namespace Array;

template<class T>
void heapsort(array1<T> arr)
{  
  // from http://www.codecodex.com/wiki/Heapsort, adapted from numerical recipes
  //GFDL license, which is compatible with LGPL, I think?
  unsigned n = arr.Size(), i = n/2, parent, child;  
  unsigned t;  
  
  for (;;) { 
    if (i > 0) {
      i--;
      t = arr[i];
    } else {
      n--;
      if (n == 0) return;
      t = arr[n];
      arr[n] = arr[0];
    }  
    
    parent = i;
    child = i*2 + 1;
    
    while (child < n) {  
      if (child + 1 < n  &&  arr[child + 1] > arr[child]) {  
	child++;
      }  
      if (arr[child] > t) {
	arr[parent] = arr[child];
	parent = child;
	child = parent*2 + 1;
      } else { 
	break;
      }  
    }  
    arr[parent] = t;
  }  
}  

template<class T>
void heapsort(DynVector<T> arr)
{  
  // from http://www.codecodex.com/wiki/Heapsort, adapted from numerical recipes
  //GFDL license, which is compatible with LGPL, I think?
  unsigned n = arr.Size(), i = n/2, parent, child;  
  unsigned t;  
  
  for (;;) { 
    if (i > 0) {
      i--;
      t = arr[i];
    } else {
      n--;
      if (n == 0) return;
      t = arr[n];
      arr[n] = arr[0];
    }  
    
    parent = i;
    child = i*2 + 1;
    
    while (child < n) {  
      if (child + 1 < n  &&  arr[child + 1] > arr[child]) {  
	child++;
      }  
      if (arr[child] > t) {
	arr[parent] = arr[child];
	parent = child;
	child = parent*2 + 1;
      } else { 
	break;
      }  
    }  
    arr[parent] = t;
  }  
}  
