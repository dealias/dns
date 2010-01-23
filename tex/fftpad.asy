int p=2,q=3;
int m=4;
int n=q*m;

pair zeta=exp(2*pi*I/n);
pair[] Zeta=new pair[n];

for(int i=0; i < n; ++i)
  Zeta[i]=zeta^i;

// return an fft with p*m data padded to n.
pair[] fftpad(pair[] f, int sign=-1)
{
  pair[] F=new pair[m];
  pair[] h=new pair[p*m];

  for(int r=0; r < q; ++r) {
    for(int k=0; k < m; ++k)
      F[k]=0.0;
    for(int a=0; a < p; ++a) {
      for(int k=0; k < m; ++k) {
        int K=k+a*m;
        F[k] += Zeta[sign*r*K % n]*f[K];
      }
    }
    F=fft(F,sign);
    for(int i=0; i < m; ++i)
      h[q*i+r]=F[i];
  }
  return h;
}

// return an fft with p*m data padded to n.
pair[] fftpadinv(pair[] f, int sign=1)
{
  pair[] F=new pair[m];
  pair[] h=new pair[p*m];

  for(int k=0; k < p*m; ++k)
    h[k]=0.0;

  for(int r=0; r < q; ++r) {
    for(int i=0; i < m; ++i)
      F[i]=f[q*i+r];

    F=fft(F,sign);
    for(int a=0; a < p; ++a) {
      for(int k=0; k < m; ++k) {
        int K=k+a*m;
        h[K] += Zeta[sign*r*K % n]*F[k];
      }
    }
  }
  return h;
}

pair[] f=sequence(p*m);
write(f);
write();

write(fftpadinv(fftpad(f,-1),1)/n);

