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
// f has length m+1, u is a work array of length m/2+1.
real[] fftpad(pair[] f, pair[] u)
{
  int m=f.length;
  int c=quotient(m,2);
  bool even=2c == m;
  assert(even);
  
  int n=3*m;

  pair zeta=exp(2*pi*I/n);
  pair[] Zeta=new pair[c+1];

  pair product=Zeta[0]=1;
  for(int i=1; i <= c; ++i) {
    product *= zeta;
    Zeta[i]=product;
  }
  
  pair zeta3c=(-0.5,-0.5*sqrt(3.0));

  real f0=f[0].x;
  u[0]=f0;
  f[m]=f0;
  for(int k=1; k <= c; ++k) {
    pair fk=f[k];
    pair fmk=conj(f[m-k]);
    u[k]=fk+fmk;
    pair zetak=Zeta[k];
    pair A=zetak*(fk.x+zeta3c*fmk.x);
    pair B=I*zetak*(fk.y+zeta3c*fmk.y);
    f[k]=A+B;
    f[m-k]=conj(A-B);
  }

  /* Unrolled loop for k=c:
  pair fk=f[c];
  f[c]=2.0*fk.x;
  pair fc=fk.x+I*fk.y*sqrt(3.0);
  u[c]=fc;
  */
  
  //  pair fc=f[c];
  real[] f2=crfft(f[c:m+1],even,-1);

  // Not necessary for convolution.
  for(int i=1; i < m; i += 2)
    f2[i]=-f2[i];

  //  f[c]=fc;
  real[] f1=crfft(f[0:c+1],even);
  real[] f0=crfft(u,even);
  real[] h=concat(f1,f2,f0);
  
  /*
  // unscramble
  real[] h=new real[3c];
  
  h[0]=f0[0];
  h[1]=f1[0];
  h[3*m-1]=f2[0];
  for(int i=1; i < m; ++i) {
    h[3*i-1]=f2[i];
    h[3*i]=f0[i];
    h[3*i+1]=f1[i];
    }*/
  
  return h;
}

// Unrolled scrambled rc version for p=2, q=3
// with n=3m; m even for now
// f has length m+1, u is a work array of length m/2+1.
pair[] fftpadinv(real[] f)
{
  int n=f.length;
  int m=quotient(n,3);
  assert(n == 3m);
  int c=quotient(m,2);
  
  pair zeta=exp(2*pi*I/n);
  pair[] Zeta=new pair[m+1];

  pair product=Zeta[0]=1;
  for(int i=1; i <= c; ++i) {
    product *= zeta;
    Zeta[i]=product;
  }
  
  pair[] f1=rcfft(f[0:m]);
  pair[] f2=rcfft(f[m:2m]);
  pair[] f0=rcfft(f[2m:n]);

  pair[] F=new pair[m];

  pair zeta3=(-0.5,0.5*sqrt(3.0));
  pair zeta3c=(-0.5,-0.5*sqrt(3.0));

  int stop=m-c-1;
  real ninv=1/n;
  F[0]=(f0[0].x+f1[0].x+f2[0].x)*ninv;
  for(int k=1; k <= stop; ++k) {
    pair f0k=f0[k];
    pair f1k=conj(Zeta[k])*f1[k];
    pair f2k=Zeta[k]*f2[k];
    F[k]=(f0k+f1k+f2k)*ninv;
    F[m-k]=conj(f0k+zeta3*f1k+zeta3c*f2k)*ninv;
  }
  bool even=true;
  
  if(even) F[c]=(f0[c].x+f1[c].x*conj(Zeta[c])+f2[c].x*Zeta[c])*ninv;

  return F;
}

pair[] convolve0(pair[] f, pair[] g)
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
  n += 2;
  int np=quotient(n,2)+1;
  
  F=copy(F);
  G=copy(G);
  
  for(int i=F.length; i < np; ++i) {
    G[i]=0.0;
    F[i]=0.0;
  }
  bool even=n % 2 == 0;
  write(n);
  
  write("F=",crfft(F,even));

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
pair[] d={-5,(3,1),(4,-2),(-3,1),(0,-2),(0,1),(4,0),(-3,-1),(1,2),(2,1),(3,1),3};
//pair[] d={-5,(3,1),(4,-2),(-3,1),(0,-2),(0,1),(4,0),(-3,-1),(1,2),(2,1)};

pair[] f=copy(d);
pair[] g=copy(d);
//pair[] g=copy(d+1);

write();

write(convolve(f,g));
write();
write(direct(f,g));
write();
write(convolve0(f,g));
write();

int c=quotient(f.length,2);
pair[] u=new pair[c+1];
//write(fftpad(f,u));
write(f);
write();

write(fftpadinv(fftpad(f,u)));
