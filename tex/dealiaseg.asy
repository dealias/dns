// general formulation of p/q padding

int p=2,q=5;
write((string) p +"/" +(string) q +" padding");
int m=4;
int n=q*m;
real eps=1e-13;
pair zeta=exp(2*pi*I/n);

pair[] f=sequence(n);

for(int i=p*m; i < n; ++i)
  f[i]=0;

write(f);
write();

pair[] g=sequence(n);

for (int r=0; r < q; ++r) {
  write();
  write("r="+(string)r);
  pair[] fm=sequence(m);
  for (int k=0; k < m; ++k) {
    for (int k=0; k < m; ++k) {
      fm[k]=0;
      for (int a=0; a < p; ++a)
	fm[k] += f[k+a*m]*zeta^(r*a*m);
      fm[k] *= zeta^(r*k);
    }
  } 
  write(fm);
  pair[] gm=fft(fm,1);
  write();
    for (int i=0; i < m; ++i) {
      g[i*q+r]=gm[i];
    }
  }

f=fft(g,-1);

// FIXME: why doesn't this work at all? asy error?
for(int i=0; i < m; ++i) {
  if(abs(f[i].x) < eps) f[i]=(0.0,f[i].y);
  if(abs(f[i].y) < eps) f[i]=f[i].x;
}
write(n);
write(f/n);

pair[] f0=sequence(n);
for(int i=p*m; i < n; ++i)
  f0[i]=0;
real error=0.0;
for(int i=0; i < n; ++i)
  error += abs(f0[i]-f[i]/n)^2;

write("error="+(string) sqrt(error));
