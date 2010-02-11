import math;
import utils;

int nx=3;
int ny=3;

real[][] f=new real[nx][ny];

for(int i=0; i < nx; ++i) 
  for(int j=0; j < ny; ++j)
    f[i][j]=i+j;
//    f[i][j]=unitrand();

write(f);
write();

real ninv=1/(nx*ny);

int n;
int m;
pair[] Zeta;

// Unrolled scrambled version for p=1, q=2.
pair[] fftpad(pair[] f)
{
  for(int k=0; k < m; ++k)
    f[k+m]=conj(Zeta[k])*f[k];
  f=concat(fft(f[0:m],-1),fft(f[m:2m],-1));
  return f;
}

// Unrolled scrambled version for p=1, q=2.
pair[] fftpadinv(pair[] f)
{
  f=concat(fft(f[0:m],1),fft(f[m:2m],1));
  for(int k=0; k < m; ++k)
    f[k] += Zeta[k]*f[k+m];
  return f[0:m];
}

pair[][] fftpad(pair[][] a)
{
  pair[][] A=new pair[a.length][];
  int k=0;
  for(pair[] v : a) {
    A[k]=fftpad(v);
    ++k;
  }
  a=transpose(A);
  k=0;
  for(pair[] v : a) {
    A[k]=fftpad(v);
    ++k;
  }
  return transpose(A);
}

pair[][] fftpadinv(pair[][] a)
{
  pair[][] A=new pair[a.length][];
  int k=0;
  for(pair[] v : a) {
    A[k]=fftpadinv(v);
    ++k;
  }
  a=transpose(A);
  k=0;
  for(pair[] v : a) {
    A[k]=fftpadinv(v);
    ++k;
  }
  return unpad(transpose(A),nx,ny);
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
pair[][] f_direct=direct(f,f);
write(f_direct);

pair[][] F=fft(pad(f,nx,ny),1);
pair[][] G=fft(pad(f,nx,ny),1);

for(int i=0; i < F.length; ++i)
  for(int j=0; j < F[0].length; ++j)
    F[i][j] *= G[i][j];

write();
real ninv=1/(F.length*F[0].length);
write("padded:");
pair[][] f_padded=ninv*unpad(fft(F,-1),nx,ny);
//write(f_padded);
write("error="+(string)maxerror(f_direct,f_padded));

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
write("fftpad(f):");
write(fftpad(f[0]));
write();


pair[][] F=fftpad(f);
pair[][] G=fftpad(f);

for(int i=0; i < F.length; ++i)
  for(int j=0; j < F[0].length; ++j)
    F[i][j] *= G[i][j];

real ninv=1/(F.length*F[0].length);
write("unpadded:");
pair[][] f_unpadded=ninv*fftpadinv(F);
//write(f_unpadded);
write("error="+(string)maxerror(f_direct,f_unpadded));
