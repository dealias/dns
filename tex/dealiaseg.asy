int p=2,q=3;
write((string) p +"/" +(string) q +" padding");
int m=4;
int n=q*m;

pair zeta=exp(-2*pi*I/n);
pair[] Zeta=new pair[n];

for(int i=0; i < n; ++i)
  Zeta[i]=zeta^i;

pair[] fp=sequence(n);

// p/q zero-padding
for(int i=p*m; i < n; ++i)
  fp[i]=0;

write("input data:");
write(fp);
pair[] gp=fft(fp,-1);
write();
write("conventional transform:");
write(gp);
fp=fft(gp,1);
write();

write("inverse transform:");
write(fp/n);

write();

write("new transform:");
pair[] f=sequence(p*m);

pair[][] g=array(q,array(m,(0,0)));
write("r="+(string) 0);
for(int a=0; a < p; ++a) {
  for(int k=0; k < m; ++k) {
    int K=k+a*m;
    g[0][k] += f[K];
  }
}
g[0]=fft(g[0],-1);
write(g[0]);

for(int r=1; r < q; ++r) {
  write("r="+(string) r);
  for(int a=0; a < p; ++a) {
    for(int k=0; k < m; ++k) {
      int K=k+a*m;
      g[r][k] += Zeta[r*K % n]*f[K];
    }
  }
  g[r]=fft(g[r],-1);
}

pair[] f=array(p*m,(0,0));
g[0]=fft(g[0],1);
for(int a=0; a < p; ++a)
  for(int k=0; k < m; ++k)
    f[k+a*m] += g[0][k];

for(int r=1; r < q; ++r) {
  g[r]=fft(g[r],1);
  for(int a=0; a < p; ++a) {
    for(int k=0; k < m; ++k) {
      int K=k+a*m;
      f[K] += Zeta[-r*K % n]*g[r][k];
    }
  }
}

write();
write(f/n);
