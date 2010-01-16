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
  int stop=c;
  if(even) --stop;
  
  int n=q*m;
  
  pair zeta=exp(-2*pi*I/n);
  pair[] Zeta=new pair[n];
  Zeta.cyclic=true;

  for(int i=0; i < n; ++i)
    Zeta[i]=zeta^i;

  write("m=",m);
  write("n=",n);

  pair[] F=new pair[m];
  pair[] G=new pair[m];

  pair[] sym(pair[] f) {
    pair[] h=new pair[c+1];
    h[0]=2f[0].x;
    for(int k=1; k <= c; ++k)
      h[k]=f[k]+conj(f[m-k]);
    return h;
  }

  // r=0:
  F=rcfft((crfft(sym(f),even)-f[0].x)*(crfft(sym(g),even)-g[0].x));

  pair[] h=new pair[m];
  h[0]=F[0].x;
  for(int k=1; k <= stop; ++k) {
    pair Fk=F[k];
    h[k]=Fk;
    h[m-k]=conj(Fk);
  }
  if(even) h[c]=F[c].x;

  // r=1:
  for(int k=0; k < m; ++k) {
    pair Zetak=Zeta[k];
    F[k]=Zetak*f[k];
    G[k]=Zetak*g[k];
  }
  
  F=rcfft((crfft(sym(F),even)-F[0].x)*(crfft(sym(G),even)-G[0].x));
    
  h[0] += F[0].x;
  pair Zetam=Zeta[m];
  for(int k=1; k <= stop; ++k) {
    pair Fk=Zeta[-k]*F[k];
    h[k] += Fk;
    h[m-k] += conj(Zetam*Fk);
  }
  if(even) h[c] += Zeta[-c]*F[c].x;

  // r=2:
  for(int k=0; k < m; ++k) {
    pair Zetamk=Zeta[-k];
    F[k]=Zetamk*f[k];
    G[k]=Zetamk*g[k];
  }

  F=rcfft((crfft(sym(F),even)-F[0].x)*(crfft(sym(G),even)-G[0].x));
    
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

//pair[] d={-5,(3,1),(4,-2),(-3,1),(0,-2),(0,1),(4,0),(-3,-1),(1,2),(2,1),(3,1)};
//pair[] d={-5,(3,1),(4,-2),(-3,1),(0,-2),(0,1),(4,0),(-3,-1),(1,2),(2,1),(3,1),3};
pair[] d={-5,(3,1),(4,-2),(-3,1),(0,-2),(0,1),(4,0),(-3,-1),(1,2),(2,1)};

pair[] f=copy(d);
pair[] g=copy(d+1);

write();

write(convolve(f,g));
write();
write(direct(f,g));
write();
write(convolve0(f,g));
