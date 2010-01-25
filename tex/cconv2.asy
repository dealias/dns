import math;

int nx=3;
int ny=3;

real[][] f=new real[nx][ny];

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

for(int i=0; i < nx; ++i) 
  for(int j=0; j < ny; ++j)
    f[i][j]=i+j;
//    f[i][j]=unitrand();

write(f);
write();

real ninv=1/(nx*ny);

pair[][] pad(pair[][] f)
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

pair[][] unpad(pair[][] f)
{
  pair[][] F=new pair[f.length][];
  for(int i=0; i < nx; ++i)
    F[i]=f[i][0:ny];
  return F;
}

int n;
int m;
pair[] Zeta;

// Unrolled scrambled version for p=1, q=2.
pair[] ffthalf(pair[] f)
{
  for(int k=0; k < m; ++k)
    f[k+m]=conj(Zeta[k])*f[k];
  f=concat(fft(f[0:m],-1),fft(f[m:2m],-1));
  return f;
}

// Unrolled scrambled version for p=1, q=2.
pair[] ffthalfinv(pair[] f)
{
  f=concat(fft(f[0:m],1),fft(f[m:2m],1));
  for(int k=0; k < m; ++k)
    f[k] += Zeta[k]*f[k+m];
  return f[0:m];
}

pair[][] ffthalf(pair[][] a)
{
  pair[][] A=new pair[a.length][];
  int k=0;
  for(pair[] v : a) {
    A[k]=ffthalf(v);
    ++k;
  }
  a=transpose(A);
  k=0;
  for(pair[] v : a) {
    A[k]=ffthalf(v);
    ++k;
  }
  return transpose(A);
}

pair[][] ffthalfinv(pair[][] a)
{
  pair[][] A=new pair[a.length][];
  int k=0;
  for(pair[] v : a) {
    A[k]=ffthalfinv(v);
    ++k;
  }
  a=transpose(A);
  k=0;
  for(pair[] v : a) {
    A[k]=ffthalfinv(v);
    ++k;
  }
  return unpad(transpose(A));
}

pair[][] direct(pair[][] F, pair[][] G)
{
  int m=F.length;
  int n=F[0].length;
  pair[][] H=new pair[m][n];
  for(int i=0; i < m; ++i) {
    for(int j=0; j < n; ++j) {
      pair sum;
      for(int k=0; k <= i; ++k)
        for(int p=0; p <= j; ++p)
        sum += F[k][p]*G[i-k][j-p];
      H[i][j]=sum;
    }
  }
  return H;
}	

write("direct:");
write(direct(f,f));

pair[][] F=fft(pad(f),1);
pair[][] G=fft(pad(f),1);

for(int i=0; i < F.length; ++i)
  for(int j=0; j < F[0].length; ++j)
    F[i][j] *= G[i][j];

write();
real ninv=1/(F.length*F[0].length);
write("padded:");
write(ninv*unpad(fft(F,-1)));

write();

assert(nx == ny);
m=nx;
n=2m;

pair zeta=exp(2*pi*I/n);
Zeta=new pair[n];

for(int i=0; i < n; ++i)
  Zeta[i]=zeta^i;


write("f:");
write(f[0]);
write();
write("ffthalf(f):");
write(ffthalf(f[0]));
write();


pair[][] F=ffthalf(f);
pair[][] G=ffthalf(f);

for(int i=0; i < F.length; ++i)
  for(int j=0; j < F[0].length; ++j)
    F[i][j] *= G[i][j];

real ninv=1/(F.length*F[0].length);
write("unpadded:");
write(ninv*ffthalfinv(F));
