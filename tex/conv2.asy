import math;
import utils;

int nx=4;
int ny=4;

real[][] f=new real[nx][ny];

real[][] crfft2d(pair[][] f, bool even=true)
{
  int nx=f.length;
  int my=f[0].length;
  int L=even ? 2my-2 : 2my-1;
  
  pair[][] h=new pair[nx][L];

  for(int i=0; i < nx; ++i)
    h[i][0]=f[i][0];
    
  for(int j=1; j < my-1; ++j) {
    h[0][j]=f[0][j];
    h[0][L-j]=conj(f[0][j]);
  }
  
  for(int i=1; i < nx; ++i) {
    for(int j=1; j < my-1; ++j) {
      h[i][j]=f[i][j];
      h[nx-i][L-j]=conj(f[i][j]);
    }
  }

  h[0][my-1]=f[0][my-1];
  if(!even) h[0][my]=conj(f[0][my-1]);

  for(int i=1; i < nx; ++i) {
    h[i][my-1]=f[i][my-1];
    if(!even) h[i][my]=conj(f[nx-i][my-1]);
  }
  
  write(h);
  write();
  h=fft(h,1);
  
  real[][] H=new real[nx][L];
  for(int i=0; i < nx; ++i)
    H[i]=map(xpart,h[i]);
  
  return H;
}

pair[][] rcfft2d(real[][] f)
{
  int nx=f.length;
  int ny=f[0].length;
  int my=quotient(ny,2)+1;
  
  pair[][] F=fft((pair[][]) f,-1);

  for(int i=0; i < F.length; ++i)
  F[i]=F[i][0:my];
  return F;
}

for(int i=0; i < nx; ++i) 
  for(int j=0; j < ny; ++j)
    f[i][j]=i+j;

write(f);
write();

real ninv=1/(nx*ny);

write(ninv*crfft2d(rcfft2d(f),ny % 2 == 0));

pair[][] pad(pair[][] f)
{
  int nx=f.length;
  int ny=f[0].length;
  
  int Nx=floor(quotient(nx*3,2));
  int Ny=floor(quotient(ny*3,2));

  pair[][] g=array(Nx,array(Ny,(0,0)));
  for(int i=0; i < nx; ++i)
    for(int j=0; j < ny; ++j)
      g[i][j]=f[i][j];
  return g;
}

pair[][] unpad(pair[][] f)
{
  pair[][] F=new pair[f.length][];
  for(int i=0; i < nx; ++i)
    F[i]=f[i][0:ny];
  return F;
}
