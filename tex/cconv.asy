int p=1,q=2;
write((string) p +"/" +(string) q +" padding");

pair[] convolve0(pair[] f, pair[] g)
{
  int m=f.length;
  int n=quotient(m*q,p);
  
  pair zeta=exp(2*pi*I/n);
  pair[] Zeta=new pair[n];

  for(int i=0; i < n; ++i)
    Zeta[i]=zeta^i;

  for(int k=0; k < m; ++k) {
    f[k+m]=conj(Zeta[k])*f[k];
    g[k+m]=conj(Zeta[k])*g[k];
  }
  
  f=concat(fft(fft(f[0:m],-1)*fft(g[0:m],-1),1),
           fft(fft(f[m:2m],-1)*fft(g[m:2m],-1),1));

  for(int k=0; k < m; ++k)
    f[k]=(f[k]+Zeta[k]*f[k+m])/n;

  return f[0:m];
}

pair[] convolve(pair[] F, pair[] G)
{
  int m=F.length;
  int n=quotient(m*q,p);
  
  F=copy(F);
  G=copy(G);
  
  for(int i=m; i < n; ++i) {
    F[i]=0.0;
    G[i]=0.0;
  }

  return (fft(fft(F,-1)*fft(G,-1),1)/n)[0:m];
}	

pair[] direct(pair[] F, pair[] G)
{
  int m=F.length;
  pair[] H=new pair[m];
  for(int i=0; i < m; ++i) {
    pair sum;
    for(int j=0; j <= i; ++j) sum += F[j]*G[i-j];
    H[i]=sum;
  }
  return H;
}	

pair[] d={(-5,3),(3,1),(4,-2),(-3,1),(0,-2),(0,1),(4,0),(-3,-1),(1,2),(2,1),(3,1)};

pair[] f=copy(d);
pair[] g=copy(d);

write();

write(convolve(f,g));
write();
write(direct(f,g));
write();
write(convolve0(f,g));
