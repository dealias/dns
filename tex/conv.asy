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

real[] crfft0(pair[] f)
{
  int m=f.length;
  int L=2m-1;
  pair[] h=new pair[L];

  h[0]=f[0];
  for(int i=1; i < m; ++i) {
    h[i]=f[i];
    h[L-i]=conj(f[i]);
  }
  return map(xpart,fft(h,-1));
}

pair[] rcfft0(real[] f)
{
  return fft(f,1)[0:quotient(f.length,2)+1];
}

// Return the non-negative Fourier components of a real vector f.
pair[] rcfft(real[] f) 
{
  return compress(fft(f,-1));
}

pair[] convolve0(pair[] f, pair[] g)
{
  int m=f.length;
  int c=quotient(m,2);
  bool even=2*c == m;
  int M;
  int start;
  if(even) {
    M=m+1;
    start=2;
  } else {
    M=m;
    start=1;
  }
  
  int n=q*M;
  
  pair zeta=exp(-2*pi*I/n);
  pair[] Zeta=new pair[n];
  Zeta.cyclic=true;

  for(int i=0; i < n; ++i)
    Zeta[i]=zeta^i;

  write("M=",M);
  write("n=",n);

  pair[] F=array(m,(0,0));
  pair[] G=array(m,(0,0));

  write("r="+(string) 0);

  pair[] sym(pair[] f) {
    pair[] h=new pair[c+1];
    h[0]=2f[0].x;
    if(even) 
      h[1]=f[1];
    for(int k=start; k <= c; ++k)
      h[k]=f[k]+conj(f[M-k]);
    return h;
  }
  
  F=rcfft0((crfft0(sym(f))-f[0].x)*(crfft0(sym(g))-g[0].x));

  pair[] h=new pair[m];
  h[0]=F[0];
  if(even)
    h[1]=F[1];
  for(int k=start; k <= c; ++k) {
    pair Fk=F[k];
    h[k]=Fk;
    h[M-k]=conj(Fk);
  }

  for(int r=1; r < q; ++r) {
    F=array(m,(0,0));
    G=array(m,(0,0));
    write("r="+(string) r);
    for(int k=0; k < m; ++k) {
      pair Zetark=Zeta[r*k];
      F[k] += Zetark*f[k];
      G[k] += Zetark*g[k];
    }

    F=rcfft0((crfft0(sym(F))-F[0].x)*(crfft0(sym(F))-F[0].x));
    h[0] += F[0];
    if(even)
      h[1] += Zeta[-r]*F[1];
    pair ZetarM=Zeta[r*M];
    for(int k=start; k <= c; ++k) {
      pair Fk=Zeta[-r*k]*F[k];
      h[k] += Fk;
      h[M-k] += conj(ZetarM*Fk);
    }
  }
  return h/n;
}

int n=64;
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
pair[] d={-5,(3,1),(4,-2),(-3,1),(0,-2),(0,1),(4,0),(-3,-1),(1,2),(2,1),(3,1),3};

pair[] f=copy(d);
pair[] g=copy(d);

write();

write(convolve(f,g));
write();
write(direct(f,g));
write();
write(convolve0(f,g));
