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

  pair[] pack(pair[] f, pair[] g) {
    pair[] h;
    h[0]=2(f[0].x,g[0].x);
    int c=quotient(m,2);
    for(int k=1; k <= c; ++k) {
      int K=m-k;
      pair F=f[k]+conj(f[K]);
      pair G=g[k]+conj(g[K]);
      h[k]=F+I*G;
      h[K]=conj(F-I*G);
    }
    return h;
  }
  
  pair[] H=fft(pack(f,g),1)-(f[0].x,g[0].x);
  pair[] h=fft(map(xpart,H)*map(ypart,H),-1)/n;

  for(int r=1; r < q; ++r) {
    F=array(m,(0,0));
    G=array(m,(0,0));
    write("r="+(string) r);
    for(int k=0; k < m; ++k) {
      F[k] += Zeta[r*k]*f[k];
      G[k] += Zeta[r*k]*g[k];
    }
    
    H=fft(pack(F,G),1)-(F[0].x,G[0].x);
    F=fft(map(xpart,H)*map(ypart,H),-1)/n;
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
