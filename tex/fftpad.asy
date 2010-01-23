int p=2,q=3;
int m=4;
int n=q*m;

pair zeta=exp(2*pi*I/n);
pair[] Zeta=new pair[n];

for(int i=0; i < n; ++i)
  Zeta[i]=zeta^i;

// return an fft with p*m data padded to n.
pair[] fftpad(pair[] f)
{
  pair[] F=new pair[m];
  pair[] h=new pair[p*m];

  for(int r=0; r < q; ++r) {
    for(int k=0; k < m; ++k)
      F[k]=0.0;
    for(int a=0; a < p; ++a) {
      for(int k=0; k < m; ++k) {
        int K=k+a*m;
        F[k] += Zeta[-r*K % n]*f[K];
      }
    }
    F=fft(F,-1);
    for(int i=0; i < m; ++i)
      h[q*i+r]=F[i];
  }
  return h;
}

// Unrolled version for p=2, q=3.
pair[] fftpad(pair[] f)
{
  pair[] H=new pair[m];
  pair[] h=new pair[p*m];

  for(int k=0; k < m; ++k) {
    pair z1=Zeta[-k % n];
    pair z2=z1*Zeta[-m % n];
    pair fk=f[k];
    pair fkm=f[k+m];
    pair a=z1.x*fk+z2.x*fkm;
    pair b=I*(z1.y*fk+z2.y*fkm);
    f[k]=a+b;
    f[k+m]=a-b;
    H[k]=fk+fkm;
  }

  f=concat(fft(f[0:m],-1),fft(f[m:2m],-1));
  H=fft(H,-1);

  h[0]=H[0];
  h[1]=f[0];
  h[n-1]=f[m];
  for(int i=1; i < m; ++i) {
    h[3*i-1]=f[m+i];
    h[3*i]=H[i];
    h[3*i+1]=f[i];
  }
  return h;
}

// return an fft with p*m data padded to n.
pair[] fftpadinv(pair[] f)
{
  pair[] F=new pair[m];
  pair[] h=new pair[p*m];

  for(int k=0; k < p*m; ++k)
    h[k]=0.0;

  for(int r=0; r < q; ++r) {
    for(int i=0; i < m; ++i)
      F[i]=f[q*i+r];

    F=fft(F,1);
    for(int a=0; a < p; ++a) {
      for(int k=0; k < m; ++k) {
        int K=k+a*m;
        h[K] += Zeta[r*K % n]*F[k];
      }
    }
  }
  return h;
}

// Unrolled scrambled version for p=2, q=3.
pair[] ffttwothirds(pair[] f)
{
  pair zetamc=conj(Zeta[m]);
  for(int k=0; k < m; ++k) {
    pair z1=conj(Zeta[k]);
    pair z2=z1*zetamc;
    pair fk=f[k];
    pair fkm=f[k+m];
    pair a=z1.x*fk+z2.x*fkm;
    pair b=I*(z1.y*fk+z2.y*fkm);
    f[k]=a+b;
    f[k+m]=a-b;
    f[k+2m]=fk+fkm;
  }

  f=concat(fft(f[0:m],-1),fft(f[m:2m],-1),concat(fft(f[2m:3m],-1)));
  return f;
}

// Unrolled scrambled version for p=2, q=3.
pair[] ffttwothirdsinv(pair[] f)
{
  f=concat(fft(f[0:m],1),fft(f[m:2m],1),concat(fft(f[2m:3m],1)));
  pair zetam=Zeta[m];
  for(int k=0; k < m; ++k) {
    pair fk=f[k];
    pair fkm=f[k+m];
    pair fkmm=f[k+2m];
    pair z1=Zeta[k];
    pair z2=z1*zetam;
    f[k]=z1*fk+conj(z1)*fkm+fkmm;
    f[k+m]=z2*fk+conj(z2)*fkm+fkmm;
  }
  return f;
}

// Unrolled scrambled version for p=1, q=2.
pair[] ffthalf(pair[] f)
{
  for(int k=0; k < m; ++k)
    f[k+m]=conj(Zeta[k])*f[k];
  f=concat(fft(f[0:m],-1),fft(f[m:2m],-1));
  return f;
}

// Unrolled scrambled version for p=2, q=3.
pair[] ffthalfinv(pair[] f)
{
  f=concat(fft(f[0:m],1),fft(f[m:2m],1));
  for(int k=0; k < m; ++k)
    f[k] += Zeta[k]*f[k+m];
  return f;
}

pair[] f=sequence(p*m);
write(f);
write();

//write((ffttwothirdsinv(ffttwothirds(f))/n)[0:p*m]);

p=1;
q=2;
m=8;
n=q*m;

pair zeta=exp(2*pi*I/n);

for(int i=0; i < n; ++i)
  Zeta[i]=zeta^i;

pair[] f=sequence(p*m);
write((ffthalfinv(ffthalf(f))/n)[0:p*m]);

