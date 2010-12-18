#include "dnsbase.h"

void DNSBase::CasimirTransfer(const vector2& Src, const vector2& Y)
{
  Set(Tn,Src[TRANSFERN]);
    
  for(unsigned K=0; K < nshells; K++)
    Tn[K]=0.0;
    
  static unsigned int my1=my+1;
  static unsigned int origin=mx*my1;

  f(origin)=0.0;
  g(origin)=0.0;
  h(origin)=0.0;
    
  for(unsigned i=0; i < Nx; ++i) {
    vector wi=w[i];
    vector fi=f[i];
    vector gi=g[i];
    vector hi=h[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      Complex wij=wi[j];
      fi[j]=wij;
      gi[j]=wij;
      hi[j]=wij;
    }
  }
    
  TConvolution->convolve(f,g,h);
    
  for(unsigned i=0; i < Nx; i++) {
    int I=(int) i-(int) xorigin;
    int I2=I*I;
    vector fi=f[i];
    vector Si=f0[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      // FIXME: optimize to work over a quadrant
      // FIXME: changes with spectrum type
      Tn[(unsigned)(sqrt(I2+j*j)-0.5)].re += realproduct(Si[j],fi[j]);
    }
  }

  for(unsigned i=0; i < Nx; ++i) {
    vector wi=w[i];
    vector fi=f[i];
    vector gi=g[i];
    vector hi=h[i];
    vector Si=f0[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      Complex wij=wi[j];
      fi[j]=Si[j];
      gi[j]=wij;
      hi[j]=wij;
    }
  }
    
  TConvolution->convolve(f,g,h);
    
  for(unsigned i=0; i < Nx; i++) {
    int I=(int) i-(int) xorigin;
    int I2=I*I;
    vector fi=f[i];
    vector wi=w[i];
    for(unsigned j=i <= xorigin ? 1 : 0; j < my; ++j) {
      // FIXME: optimize to work over a quadrant
      // FIXME: changes with spectrum type
      Tn[(unsigned)(sqrt(I2+j*j)-0.5)].re += 3.0*realproduct(wi[j],fi[j]);
    }
  }
}
