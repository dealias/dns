int p=2,q=3;
write((string) p +"/" +(string) q +" padding");

// Return the inverse Fourier transform of a Hermitian vector f.
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

  pair[] F=array(m,(0,0));
  pair[] G=array(m,(0,0));

  pair[] sym(pair[] f) {
    pair[] h=new pair[c+1];
    h[0]=2f[0].x;
    for(int k=1; k <= c; ++k)
      h[k]=f[k]+conj(f[m-k]);
    return h;
  }

  F=rcfft((crfft(sym(f),even)-f[0].x)*(crfft(sym(g),even)-g[0].x));

  pair[] h=new pair[m];
  h[0]=F[0].x;
  for(int k=1; k <= stop; ++k) {
    pair Fk=F[k];
    h[k]=Fk;
    h[m-k]=conj(Fk);
  }
  if(even) h[c]=F[c].x;

  for(int r=1; r < q; ++r) {
    F=array(m,(0,0));
    G=array(m,(0,0));
    for(int k=0; k < m; ++k) {
      pair Zetark=Zeta[r*k];
      F[k] += Zetark*f[k];
      G[k] += Zetark*g[k];
    }

    F=rcfft((crfft(sym(F),even)-F[0].x)*(crfft(sym(F),even)-F[0].x));
    
    h[0] += F[0].x;
    pair Zetarm=Zeta[r*m];
    for(int k=1; k <= stop; ++k) {
      pair Fk=Zeta[-r*k]*F[k];
      h[k] += Fk;
      h[m-k] += conj(Zetarm*Fk);
    }
    if(even) h[c] += Zeta[-r*c]*F[c].x;
  }
  return h/n;
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
  
  return rcfft((crfft(F)*crfft(G))/n)[0:m];
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
pair[] g=copy(d);

write();

write(convolve(f,g));
write();
write(direct(f,g));
write();
write(convolve0(f,g));
