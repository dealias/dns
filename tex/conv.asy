int p=2,q=3;
write((string) p +"/" +(string) q +" padding");
int m=4;
int n=q*m;

pair zeta=exp(-2*pi*I/n);
pair[] Zeta=new pair[n];

for(int i=0; i < n; ++i)
  Zeta[i]=zeta^i;

//pair d[]={-5,(3,1),(4,-2),(-3,1),(0,-2),(0,1),(4,0),(-3,-1),(1,2),(2,1),(3,1)};
pair[] d=sequence(p*m);

pair[] f=copy(d);
pair[] g=copy(d);

pair[] F=array(m,(0,0));
pair[] G=array(m,(0,0));

pair[] h=array(p*m,(0,0));

pair d[];
for(int i=0; i < p*m; ++i) {
  pair sum;
  for(int j=0; j <= i; ++j) sum += f[j]*g[i-j];
  for(int j=i+1; j < p*m; ++j) sum += f[j]*conj(g[j-i]);
  for(int j=1; j < p*m-i; ++j) sum += conj(f[j])*g[i+j];
  d[i]=sum;
}	

write(d);


void shift(pair[] f) 
{
  for(int i=1; i < f.length; i += 2)
    f[i]=-f[i];
}

write("r="+(string) 0);

shift(f);
shift(g);

for(int a=0; a < p; ++a) {
  for(int k=0; k < m; ++k) {
    int K=k+a*m;
    F[k] += f[K];
    G[k] += g[K];
  }
}

F=fft(fft(F,-1)*fft(G,-1),1)/n;

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
  F=fft(fft(F,-1)*fft(G,-1),1)/n;
  for(int a=0; a < p; ++a) {
    for(int k=0; k < m; ++k) {
      int K=k+a*m;
      h[K] += conj(Zeta[r*K % n])*F[k];
    }
  }
}

shift(h);

write();
write(h);

