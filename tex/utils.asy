pair[][] operator *(pair[][] a, pair b)
{
  int n=a.length;
  pair[][] m=new pair[n][];
  for(int i=0; i < n; ++i)
    m[i]=a[i]*b;
  return m;

}

pair[][] operator *(pair b, pair[][] a)
{
  return a*b;
}

pair[][] operator cast(real[][] f) 
{
  pair[][] F=new pair[f.length][];
  for(int i=0; i < f.length; ++i)
    F[i]=(pair[]) f[i];
  return F;
}

real maxerror(pair[][] f, pair[][] g) {
  real e=0;
  for (int i=0; i < f.length; ++i) {
    for (int j=0; j < f[i].length; ++j) {
      if (e < abs(f[i][j]-g[i][j]))
	e=abs(f[i][j]-g[i][j]);
    }
  }
  return e;
}

pair[][] pad(pair[][] f,int nx, int ny)
{
  int nx=f.length;
  int ny=f[0].length;
  
  int Nx=2*nx;
  int Ny=2*ny;
  //  int Nx=floor(quotient(nx*3,2));
  //  int Ny=floor(quotient(ny*3,2));

  pair[][] g=array(Nx,array(Ny,(0,0)));
  for(int i=0; i < nx; ++i)
    for(int j=0; j < ny; ++j)
      g[i][j]=f[i][j];
  return g;
}

pair[][] unpad(pair[][] f,int nx, int ny)
{
  pair[][] F=new pair[f.length][];
  for(int i=0; i < nx; ++i)
    F[i]=f[i][0:ny];
  return F;
}
