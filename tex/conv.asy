int p=2,q=3;
write((string) p +"/" +(string) q +" padding");
int m=4;
int n=q*m;

pair zeta=exp(-2*pi*I/n);
pair[] Zeta=new pair[n];

for(int i=0; i < n; ++i)
  Zeta[i]=zeta^i;

//pair d[]={-5,(3,1),(4,-2),(-3,1),(0,-2),(0,1),(4,0),(-3,-1),(1,2),(2,1),(3,1)};

int pm=p*m;

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

// Return the inverse Fourier transform of a Hermitian vector f.
real[] crfft(pair[] f) 
{
  return map(xpart,shift(fft(concat(map(conj,reverse(f)[0:f.length-1]),
                                    f[0:f.length-1]),1)));
}

// Return the non-negative Fourier components of a real vector f.
pair[] rcfft(real[] f) 
{
  return fft(f,-1)[quotient(f.length,2):f.length];
}

/*
pair[] e=new pair[] {(1,2),(3,4),(5,6),(7,8)};
pair[] D=new pair[2*m];

D[0]=0;
for(int i=1; i < m; ++i)
  D[i]=conj(e[m-i]);
D[m]=7;
for(int i=1; i < m; ++i)
  D[m+i]=e[i];

write(D);

write();

//d[pm-1]=0;

pair[] d=D[4:8];
write(d);

write();

write(crfft(d));
write();
write(rcfft(crfft(d))/pm);
write();

write(fft(D,-1));



pair[] f=copy(d);
pair[] g=copy(d);

pair[] F=array(m,(0,0));
pair[] G=array(m,(0,0));

pair[] h=array(pm,(0,0));
*/

/*
pair d[];
for(int i=0; i < pm; ++i) {
  pair sum;
  for(int j=max(i-3,0); j <= min(i+3,6); ++j) sum += f[j]*g[3+i-j];
  //  for(int j=0; j <= i; ++j) sum += f[j]*g[i-j];
  d[i]=sum;
}
write(d);
*/

/*
int N=n;
pair[] fp=array(N,(0,0));
pair[] gp=array(N,(0,0));
for(int i=0; i < pm; ++i) {
  fp[i]=f[i];
  gp[i]=g[i];
}

write();
pair[] P=fft(fp,-1)*fft(gp,-1);
P=shift(P);
write((fft(P,1)/N)[0:pm]);
*/

/*
write("r="+(string) 0);

for(int a=0; a < p; ++a) {
  for(int k=0; k < m; ++k) {
    int K=k+a*m;
    F[k] += f[K];
    G[k] += g[K];
  }
}

F=fft(shift(fft(F,-1)*fft(G,-1)),1)/n;

for(int a=0; a < p; ++a)
  for(int k=0; k < m; ++k)
    h[k+a*m] += F[k];

for(int r=1; r < q; ++r) {
  F=array(m,(0,0));
  G=array(m,(0,0));
  write("r="+(string) r);
  for(int a=0; a < p; ++a) {
    for(int k=0; k < m; ++k) {
      int K=k+a*m;
      F[k] += Zeta[r*K % n]*f[K];
      G[k] += Zeta[r*K % n]*g[K];
    }
  }
  F=fft(shift(fft(F,-1)*fft(G,-1)),1)/n;
  for(int a=0; a < p; ++a) {
    for(int k=0; k < m; ++k) {
      int K=k+a*m;
      h[K] += conj(Zeta[r*K % n])*F[k];
    }
  }
}

write();
write(h);
*/

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
