int p=2,q=3;
write((string) p +"/" +(string) q +" padding");

// Return the inverse Fourier transform of size n for a Hermitian vector f
// of length (n/2+1). The flag even indicates whether n is even.
real[] crfft(pair[] f, bool even=true, int sign=1)
{
  int m=f.length;
  int L=even ? 2m-2 : 2m-1;
  pair[] h=new pair[L];
  h[0]=f[0];
  for(int i=1; i < m-1; ++i) {
    h[i]=f[i];
    h[L-i]=conj(f[i]);
  }
  h[m-1]=f[m-1];
  if(!even) h[m]=conj(f[m-1]);
  return map(xpart,fft(h,sign));
}

// Return the non-negative Fourier spectrum of a real vector f.
pair[] rcfft(real[] f, int sign=-1)
{
  return fft(f,sign)[0:quotient(f.length,2)+1];
}

// Unrolled scrambled cr version for p=2, q=3
// with n=3m; m even for now
// f has length m, u is a work array of length m/2+1.
real[] fftpad(pair[] f, pair[] u, bool unscramble=true)
{
  int m=f.length;
  int c=quotient(m,2);
  bool even=2c == m;
  assert(even);
  
  int n=3*m;

  pair zeta=exp(2*pi*I/n);
  pair zeta3c=(-0.5,-0.5*sqrt(3.0));

  real f0=f[0].x;
  u[0]=f0;
  pair fc=f[c];
  int m1=m-1;
  pair fmk=conj(f[m1]);
  f[m1]=f0;
  pair Zetak=zeta;
  for(int k=1; k < c; ++k) {
    pair fk=f[k];
    f[k]=fk+fmk;
    pair A=Zetak*(fk.x+zeta3c*fmk.x);
    pair B=I*Zetak*(fk.y+zeta3c*fmk.y);
    Zetak *= zeta;
    u[k]=conj(A-B);
    int mk=m1-k;
    fmk=conj(f[mk]);
    f[mk]=A+B;
  }
  
  real A=fc.x;
  real B=sqrt(3)*fc.y;
  fc=f[c];
  f[c]=2.0*A;
  u[c]=A+B;
  A -= B;

  real[] f0=crfft(f[0:c+1],even);
  //  real fcm1=f[c-1];

  f[c-1]=A;
  f[c]=fc;
  real[] f1=crfft(f[c-1:m],even,-1);
  //  f[m-1]=fcm1; // This is where we will store f0[c-1] in scrambled format.

  // Not necessary for convolution.
  for(int i=1; i < m; i += 2)
    f1[i]=-f1[i];

  real[] f2=crfft(u,even);
  
  real[] h;
  
  if(unscramble) {
    h=new real[3c];
  
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

// Unrolled scrambled rc version for p=2, q=3
// with n=3m
// f has length 3m.
pair[] fftpadinv(real[] f, bool unscramble=true)
{
  assert(!unscramble); // Not yet implemented
  
  int n=f.length;
  int m=quotient(n,3);
  assert(n == 3m);
  int c=quotient(m,2);
  bool even=2c == m;
  
  pair zeta=exp(2*pi*I/n);
  
  pair[] f0=rcfft(f[0:m]);
  pair[] f1=rcfft(f[m:2m]);
  pair[] f2=rcfft(f[2m:n]);

  pair[] F=new pair[m];

  pair zeta3=(-0.5,0.5*sqrt(3.0));
  pair zeta3c=conj(zeta3);

  int stop=m-c-1;
  real ninv=1/n;
  F[0]=(f0[0].x+f1[0].x+f2[0].x)*ninv;
  pair Zetak=zeta*ninv;
  for(int k=1; k <= stop; ++k) {
    pair f0k=f0[k]*ninv;
    pair f1k=conj(Zetak)*f1[k];
    pair f2k=Zetak*f2[k];
    Zetak *= zeta;
    F[k]=f0k+f1k+f2k;
    F[m-k]=conj(f0k)+zeta3c*conj(f1k)+zeta3*conj(f2k);
  }
  
  if(even) F[c]=(f0[c].x-f1[c].x*zeta3-f2[c].x*conj(zeta3))*ninv;

  return F;
}

// f and g have length m.
// u and v are work arrays each of length m/2+1.
pair[] convolve0(pair[] f, pair[] g, pair[] u, pair[] v)
{
  int m=f.length;
  int c=quotient(m,2);
  bool even=2c == m;
  assert(even);

  int n=3*m;
  
  pair zeta=exp(2*pi*I/n);
  pair zetac=conj(zeta);
  pair zeta3=(-0.5,0.5*sqrt(3.0));

  real f0=f[0].x;
  real g0=g[0].x;
  u[0]=f0;
  v[0]=g0;
  pair fc=f[c];
  int m1=m-1;
  pair fmk=conj(f[m1]);
  f[m1]=f0;
  pair gc=g[c];
  pair gmk=conj(g[m1]);
  g[m1]=g0;
  pair Zetak=zetac;
  for(int k=1; k < c; ++k) {
    pair fk=f[k];
    f[k]=fk+fmk;
    pair A=Zetak*(fk.x+zeta3*fmk.x);
    pair B=-I*Zetak*(fk.y+zeta3*fmk.y);
    u[k]=A-B;
    int mk=m1-k;
    fmk=conj(f[mk]);
    f[mk]=A+B; // Store conjugate of desired quantity in reverse order.

    pair gk=g[k];
    g[k]=gk+gmk;
    A=Zetak*(gk.x+zeta3*gmk.x);
    B=-I*Zetak*(gk.y+zeta3*gmk.y);
    Zetak *= zetac;
    v[k]=A-B;
    gmk=conj(g[mk]);
    g[mk]=A+B;
  }
  

  real A=fc.x;
  real B=sqrt(3)*fc.y;
  fc=f[c];
  f[c]=2.0*A;
  u[c]=A+B;
  A -= B;

  real[] f0=crfft(f[0:c+1],even);

  real C=gc.x;
  B=sqrt(3)*gc.y;
  gc=g[c];
  g[c]=2.0*C;
  v[c]=C+B;
  C -= B;
  real[] g0=crfft(g[0:c+1],even);
  f0 *= g0;
  pair[] f0=rcfft(f0);
  pair overlap0=f0[c-1];
  real overlap1=f0[c].x;

  f[c-1]=A;
  f[c]=fc;
  real[] f1=crfft(f[c-1:m],even);
  g[c-1]=C;
  g[c]=gc;
  real[] g1=crfft(g[c-1:m],even);
  f1 *= g1;
  pair[] f1=rcfft(f1);
  // Data is shifted down by 1 complex.
  

  real[] f2=crfft(u,even);
  real[] g2=crfft(v,even);
  f2 *= g2;
  pair[] f2=rcfft(f2);

  pair[] F=new pair[m];

  int stop=m-c-1;
  real ninv=1/n;
  F[0]=(f0[0].x+f1[0].x+f2[0].x)*ninv;
  pair zeta3c=conj(zeta3);
  pair Zetak=zeta*ninv;

  for(int k=1; k < stop; ++k) {
    pair f0k=f0[k]*ninv;
    pair f1k=conj(Zetak)*f1[k];
    pair f2k=Zetak*f2[k];
    Zetak *= zeta;
    F[k]=f0k+f1k+f2k;
    F[m-k]=conj(f0k)+zeta3c*conj(f1k)+zeta3*conj(f2k);
  }

  assert(even);
  
  pair f0k=overlap0*ninv;
  pair f1k=conj(Zetak)*f1[c-1];
  pair f2k=Zetak*f2[c-1];
  F[c-1]=f0k+f1k+f2k;
  F[c+1]=conj(f0k)+zeta3c*conj(f1k)+zeta3*conj(f2k);
  
  //  if(even) F[c]=(f0[c].x-f1[c].x*zeta3-f2[c].x*conj(zeta3))*ninv;
  if(even) F[c]=(overlap1-f1[c].x*zeta3-f2[c].x*conj(zeta3))*ninv;

  return F;
}

// f and g have length m.
pair[] convolveold(pair[] f, pair[] g)
{
  int m=f.length;
  int c=quotient(m,2);
  bool even=2c == m;
  int n=q*m;
  
  pair zeta=exp(2*pi*I/n);
  pair[] Zeta=new pair[c+1];

  pair product=Zeta[0]=1;
  for(int i=1; i <= c; ++i) {
    product *= zeta;
    Zeta[i]=product;
  }
  
  pair Zetam=(-0.5,0.5*sqrt(3.0));
  pair Zetamc=conj(Zetam);

  write("m=",m);
  write("n=",n);

  pair[] F=new pair[c+1];
  pair[] G=new pair[c+1];
  pair[] B=new pair[c+1];

  // r=0:
  F[0]=f[0].x;
  G[0]=g[0].x;
  for(int k=1; k <= c; ++k) {
    F[k]=f[k]+conj(f[m-k]);
    G[k]=g[k]+conj(g[m-k]);
  }
  
  pair[] h=rcfft(crfft(F,even)*crfft(G,even));

  // r=1:
  for(int k=1; k <= c; ++k) {
    pair Zetak=Zeta[k];
    F[k]=Zetak*(f[k]+conj(Zetam*f[m-k]));
    G[k]=Zetak*(g[k]+conj(Zetam*g[m-k]));
  }
  
  B=rcfft(crfft(F,even)*crfft(G,even));
    
  // r=-1:
  for(int k=1; k <= c; ++k) {
    pair Zetakc=conj(Zeta[k]);
    F[k]=Zetakc*(f[k]+Zetam*conj(f[m-k]));
    G[k]=Zetakc*(g[k]+Zetam*conj(g[m-k]));
  }

  F=rcfft(crfft(F,even)*crfft(G,even));
    
  int stop=m-c-1;
  real ninv=1/n;
  h[0]=(h[0].x+B[0].x+F[0].x)*ninv;
  for(int k=1; k <= stop; ++k) {
    pair Bk=conj(Zeta[k])*B[k];
    pair Fk=Zeta[k]*F[k];
    pair hk=h[k];
    h[m-k]=conj(hk+Zetam*Bk+Zetamc*Fk)*ninv;
    h[k]=(hk+Bk+Fk)*ninv;
  }
  if(even) h[c]=(h[c].x+B[c].x*conj(Zeta[c])+F[c].x*Zeta[c])*ninv;
  return h;
}

pair[] convolve(pair[] F, pair[] G)
{
  int m=F.length;
  int n=quotient((2*m-1)*q,p);
  //  n += 2;
  int np=quotient(n,2)+1;
  
  F=copy(F);
  G=copy(G);
  
  for(int i=F.length; i < np; ++i) {
    G[i]=0.0;
    F[i]=0.0;
  }
  bool even=n % 2 == 0;
  //  write(n);
  
  //  write("F=",crfft(F,even));

  return rcfft((crfft(F,even)*crfft(G,even))/n)[0:m];
}	

pair[] direct(pair[] F, pair[] G)
{
  int m=F.length;
  pair[] H=new pair[m];
  for(int i=0; i < m; ++i) {
    pair sum;
    for(int j=0; j <= i; ++j) sum += F[j]*G[i-j];
    for(int j=i+1; j < m; ++j) sum += F[j]*conj(G[j-i]);
    for(int j=1; j < m-i; ++j) sum += conj(F[j])*G[i+j];
    H[i]=sum;
  }
  return H;
}	

//pair[] d={-5,(3,1),(4,-2),(-3,1),(0,-2),(0,1),(4,0),(-3,-1),(1,2),(2,1),(3,1)};
pair[] d={-5,(3,1),(4,-2),(-3,1),(0,-2),(0,1),(4,1),(-3,-1),(1,2),(2,1),(3,1),3};
//pair[] d={-5,(3,1),(4,-2),(-3,1),(0,-2),(0,1),(4,0),(-3,-1),(1,2),(2,1)};

pair[] f=copy(d);
pair[] g=copy(d);
//pair[] g=copy(d+1);

write();

int c=quotient(f.length,2);

pair[] u=new pair[c+1];
pair[] v=new pair[c+1];

write(convolve(f,g));
write();
write(direct(f,g));
write();
write(convolve0(f,g,u,v));
write();

write();
//write(f);
write();
//write(fftpad(f,u));
//write(fftpadinv(fftpad(f,u,false),false));
