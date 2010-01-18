int p=2,q=3;
write((string) p +"/" +(string) q +" padding");

// Return the inverse Fourier transform of size n for a Hermitian vector f
// of length (n/2+1). The flag even indicates whether n is even.
real[] crfft(pair[] f, bool even=true)
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
  return map(xpart,fft(h,1));
}

// Return the non-negative Fourier spectrum of a real vector f.
pair[] rcfft(real[] f)
{
  return fft(f,-1)[0:quotient(f.length,2)+1];
}

pair[] convolve0(pair[] f, pair[] g)
{
  int m=f.length;
  int c=quotient(m,2);
  
  bool even=2c == m;
  int stop=even ? c-1 : c;
  
  int n=q*m;
  
  pair zeta=exp(2*pi*I/n);
  pair[] Zeta=new pair[c+1];

  pair product=Zeta[0]=1;
  for(int i=1; i <= c; ++i) {
    product *= zeta;
    Zeta[i]=product;
  }
  
  pair Zetam=(-0.5,0.5*sqrt(3.0));

  write("m=",m);
  write("n=",n);

  real f0=f[0].x;
  real g0=g[0].x;

  pair[] F=new pair[c+1];
  pair[] G=new pair[c+1];

  F[0]=2f0;
  G[0]=2g0;
  for(int k=1; k <= c; ++k) {
    F[k]=f[k]+conj(f[m-k]);
    G[k]=g[k]+conj(g[m-k]);
  }
  
  // r=0:
  pair[] h=rcfft((crfft(F,even)-f0)*(crfft(G,even)-g0));

  for(int k=1; k <= stop; ++k)
    h[m-k]=conj(h[k]);

  // r=1:
  F[0]=2f0;
  G[0]=2g0;
  for(int k=1; k <= c; ++k) {
    pair Zetak=Zeta[k];
    F[k]=Zetak*(f[k]+conj(Zetam*f[m-k]));
    G[k]=Zetak*(g[k]+conj(Zetam*g[m-k]));
  }
  
  F=rcfft((crfft(F,even)-f0)*(crfft(G,even)-g0));
    
  h[0] += F[0].x;
  for(int k=1; k <= stop; ++k) {
    pair Fk=conj(Zeta[k])*F[k];
    h[k] += Fk;
    h[m-k] += conj(Zetam*Fk);
  }
  if(even) h[c] += conj(Zeta[c])*F[c].x;

  // r=2:
  F[0]=2f0;
  G[0]=2g0;
  for(int k=1; k <= c; ++k) {
    pair Zetamk=conj(Zeta[k]);
    F[k]=Zetamk*(f[k]+Zetam*conj(f[m-k]));
    G[k]=Zetamk*(g[k]+Zetam*conj(g[m-k]));
  }

  F=rcfft((crfft(F,even)-f0)*(crfft(G,even)-g0));
    
  h[0] += F[0].x;
  for(int k=1; k <= stop; ++k) {
    pair Fk=Zeta[k]*F[k];
    h[k] += Fk;
    h[m-k] += Zetam*conj(Fk);
  }
  if(even) h[c] += Zeta[c]*F[c].x;

  return h/n;
}

pair[] convolve(pair[] F, pair[] G)
{
  int m=F.length;
  int n=quotient((2*m-1)*q,p);
  int np=quotient(n,2)+1;
  
  F=copy(F);
  G=copy(G);
  
  for(int i=F.length; i < np; ++i) {
    G[i]=0.0;
    F[i]=0.0;
  }
  bool even=n % 2 == 0;
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

pair[] d={-5,(3,1),(4,-2),(-3,1),(0,-2),(0,1),(4,0),(-3,-1),(1,2),(2,1),(3,1)};
//pair[] d={-5,(3,1),(4,-2),(-3,1),(0,-2),(0,1),(4,0),(-3,-1),(1,2),(2,1),(3,1),3};
//pair[] d={-5,(3,1),(4,-2),(-3,1),(0,-2),(0,1),(4,0),(-3,-1),(1,2),(2,1)};

pair[] f=copy(d);
//pair[] g=copy(d+1);
pair[] g=copy(d);

write();

write(convolve(f,g));
write();
write(direct(f,g));
write();
write(convolve0(f,g));
