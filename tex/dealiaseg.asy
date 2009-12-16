int p=2,q=3;
write((string) p +"/" +(string) q +" padding");
int m=4;
int n=q*m;

pair zeta=exp(-2*pi*I/n);
pair[][] Zeta=array(q,array(p*m,(0,0)));

for(int r=0; r < q; ++r)
    for(int k=0; k < p*m; ++k)
      Zeta[r][k]=zeta^(r*k);

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
for(int r=0; r < q; ++r) {
  write("r="+(string) r);
  for(int a=0; a < p; ++a) {
    for(int k=0; k < m; ++k) {
      int K=k+a*m;
      g[r][k] += Zeta[r][K]*f[K];
    }
  }
  g[r]=fft(g[r],-1);
  write(g[r]);
}

pair[] conj(pair[] z) {return map(conj,z);} // TODO: Move to asy.

pair[] f=array(p*m,(0,0));
for(int r=0; r < q; ++r) {
  g[r]=fft(g[r],1);
  for(int a=0; a < p; ++a) {
    for(int k=0; k < m; ++k) {
      int K=k+a*m;
      f[K] += conj(Zeta[r][K])*g[r][k];
    }
  }
}

write();
write(f/n);
