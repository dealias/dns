pair[] convolve0(pair[] f, pair[] g)
{
  int p=1,q=2;

  int m=f.length;
  int n=quotient(m*q,p);
  
  pair zeta=exp(2*pi*I/n);

  pair Zetak=1.0;
  for(int k=0; k < m; ++k) {
    f[k+m]=conj(Zetak)*f[k];
    g[k+m]=conj(Zetak)*g[k];
    Zetak *= zeta;
  }
  
  f=concat(fft(fft(f[0:m],-1)*fft(g[0:m],-1),1),
           fft(fft(f[m:2m],-1)*fft(g[m:2m],-1),1));

  Zetak=1.0;
  for(int k=0; k < m; ++k) {
    f[k]=(f[k]+Zetak*f[k+m])/n;
    Zetak *= zeta;
  }

  return f[0:m];
}

// Unrolled scrambled shifted cc version for p=2, q=3, with n=3m.
// f has length 2m-1 with the origin at index m
// (i.e. physical wavenumber k=-m+1 to k=m-1).
// u is a work array of length m.
pair[] cfftpad(pair[] f, pair[] u, bool unscramble=true)
{
  int m=quotient(f.length+1,2);
  assert(2m-1 == f.length);
  int n=3*m;
  bool even=m % 2 == 0;

  pair zeta=exp(2*pi*I/n);
  pair zeta3c=(-0.5,-0.5*sqrt(3.0));

  pair Zetak=zeta;
  pair fk=f[0];
  int m1=m-1;
  pair fm1=f[m1];
  u[0]=f[0]=fm1;
  for(int k=1; k < m-1; ++k) {
    int mk=m1+k;
    pair fmk=f[mk];
    pair C=fk+fmk;
    pair A=Zetak*(fmk.x+zeta3c*fk.x);
    pair B=I*Zetak*(fmk.y+zeta3c*fk.y);
    Zetak *= zeta;
    fk=f[k];
    f[k]=C;
    f[mk]=A+B;
    u[k]=conj(A-B);
  }

  int k=m1;
  int mk=m1+k;
  pair fmk=f[mk];
  pair C=fk+fmk;
  pair A=Zetak*(fmk.x+zeta3c*fk.x);
  pair B=I*Zetak*(fmk.y+zeta3c*fk.y);
  f[k]=C;
  f[mk]=A+B;
  u[k]=conj(A-B);
  
  pair[] f0=fft(f[0:m]);
  f[m-1]=fm1;
  
  pair[] f1=fft(f[m-1:2m-1]);
  pair[] f2=fft(u);
  
  pair[] h;
  
  if(unscramble) {
    h=new pair[n];
  
    h[0]=f0[0];
    h[1]=f1[0];
    h[3*m-1]=f2[0];
    for(int i=1; i < m; ++i) {
      h[3*i-1]=f2[i];
      h[3*i]=f0[i];
      h[3*i+1]=f1[i];
    }
  } else h=concat(f0,f1,f2);
  
  return h;
}

pair[] convolve(pair[] F, pair[] G)
{
  int p=1,q=2;

  int m=F.length;
  int n=quotient(m*q,p);
  
  F=copy(F);
  G=copy(G);
  
  for(int i=m; i < n; ++i) {
    F[i]=0.0;
    G[i]=0.0;
  }

  return (fft(fft(F,-1)*fft(G,-1),1)/n)[0:m];
}	

pair[] direct(pair[] F, pair[] G)
{
  int m=F.length;
  pair[] H=new pair[m];
  for(int i=0; i < m; ++i) {
    pair sum;
    for(int j=0; j <= i; ++j) sum += F[j]*G[i-j];
    H[i]=sum;
  }
  return H;
}	

pair[] d={(-5,3),(3,1),(4,-2),(-3,1),(0,-2),(0,1),(4,0),(-3,-1),(1,2),(2,1),(3,1)};

pair[] f=copy(d);
pair[] g=copy(d);

write();

write(convolve(f,g));
write();
write(direct(f,g));
write();
write(convolve0(f,g));

write();
f=copy(d);
write(f);
write();

int m=quotient(f.length+1,2);
assert(2m-1 == f.length);

int c=quotient(m,2);
pair[] F=concat(f[m-1:2m-1],array(7,(0,0)),f[0:m-1]);
//write(F);

F=fft(F);

write(F);
write();


int m=quotient(f.length+1,2);
pair[] u=new pair[m];

int p=2,q=3;
write(cfftpad(f,u));

