int p=2,q=3;
write((string) p +"/" +(string) q +" padding");

private pair[] shift(pair[] f) 
{
  for(int i=1; i < f.length; i += 2)
    f[i]=-f[i];
  return f;
}

private real[] shift(real[] f)
{
  for(int i=1; i < f.length; i += 2)
    f[i]=-f[i];
  return f;
}

pair[] decompress(pair[] f)
{
  return concat(map(conj,reverse(f)[0:f.length-1]),f[0:f.length-1]);
}

pair[] compress(pair[] f)
{
  return f[quotient(f.length,2):f.length];
}

// Return the inverse Fourier transform of a Hermitian vector f.
real[] crfft(pair[] f)
{
  return map(xpart,shift(fft(decompress(f),1)));
}

// Return the non-negative Fourier components of a real vector f.
pair[] rcfft(real[] f) 
{
  return compress(fft(f,-1));
}

pair[] convolve0(pair[] f, pair[] g)
{
  int m=f.length;
  int n=q*m;
  
  pair zeta=exp(-2*pi*I/n);
  pair[] Zeta=new pair[n];
  Zeta.cyclic=true;

  for(int i=0; i < n; ++i)
    Zeta[i]=zeta^i;

  write("m=",m);
  write("n=",n);

  pair[] F=array(m,(0,0));
  pair[] G=array(m,(0,0));

  write("r="+(string) 0);
  
  pair[] F=sequence(new pair(int k) {return f[k]+conj(f[(m-k) % m]);},f.length);
  pair[] G=sequence(new pair(int k) {return g[k]+conj(g[(m-k) % m]);},f.length);
  pair[] H=fft(F+I*G,1);
  real[] Fr=map(xpart,H)-xpart(f[0]);
  real[] Gr=map(ypart,H)-xpart(g[0]);
  pair[] h=fft(Fr*Gr,-1)/n;

  for(int r=1; r < q; ++r) {
    F=array(m,(0,0));
    G=array(m,(0,0));
    write("r="+(string) r);
    for(int k=0; k < m; ++k) {
      F[k] += Zeta[r*k]*f[k];
      G[k] += Zeta[r*k]*g[k];
    }
    
    real[] Fr=2*map(xpart,fft(F,1))-xpart(F[0]);
    real[] Gr=2*map(xpart,fft(G,1))-xpart(G[0]);
    F=fft(Fr*Gr,-1)/n;
    for(int k=0; k < m; ++k)
      h[k] += Zeta[-r*k]*F[k];
  }
  return h;
}

int n=32;
int np=quotient(n,2)+1;

pair[] convolve(pair[] F, pair[] G)
{
  int m=F.length;
  
  F=copy(F);
  G=copy(G);
  
  for(int i=F.length; i < np; ++i) {
    G[i]=0.0;
    F[i]=0.0;
  }
  
  return rcfft(shift(crfft(F)*crfft(G))/n)[0:m];
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

pair[] f=copy(d);
pair[] g=copy(d);

write();

write(convolve(f,g));
write();
write(direct(f,g));
write();
write(convolve0(f,g));
