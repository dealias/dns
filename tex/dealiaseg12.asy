int m=4, n=2*m;

pair zeta=exp(-2*pi*I/n);
pair[] zetai;
for(int i=0; i < m; ++i) zetai[i]=zeta^i;

pair[] fp=sequence(n);
// "2/4" zero-padding
for(int i=m; i < n; ++i)
  fp[i]=0;
write("input data:");
write(fp);
pair[] gp=fft(fp,-1);
write();
write("conventional transform:");
write(gp);
fp=fft(gp,1);
write();
real eps=1e-14;
// get rid of numerical error
fp=abs(fp) < eps ? array(n,0) : fp;
write("inverse transform:");
write(fp/n);

write();
write("new transform:");
pair[] f=sequence(m);

write("even terms:");
pair[] ge=fft(f,-1);
write(ge);

write();
write("odd terms:");
pair[] go=fft(f*zetai,-1);
write(go);

pair[] conj(pair[] z) {return map(conj,z);} // TODO: Move to asy.

pair[] f=fft(ge,1)+conj(zetai)*fft(go,1);

// get rid of numerical error
for(int i=0; i < m; ++i) {
  if(abs(f[i].x) < eps) f[i]=(0.0,f[i].y);
  if(abs(f[i].y) < eps) f[i]=f[i].x;
}

write();
write(f/n);
