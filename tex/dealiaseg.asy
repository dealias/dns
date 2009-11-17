int m=4, n=2*m;

write("conventional transform:");
pair[] fp=sequence(n);
// "2/4" zero-padding
for(int i=m; i < n; ++i)
  fp[i]=0;
write(fp);
pair[] gp=fft(fp,-1);
write();
write(gp);
fp=fft(gp,1);
write();
// get rid of numerical error
fp=abs(fp) < 1e-15 ? array(n,0) : fp;
write(fp/n);

write();
write("new transform:");
pair[] f=sequence(m);

write("even terms:");
write(f);
pair[] ge=fft(f,-1);
pair zetan=exp(-2*pi*I/n);
for(int i=0; i < m; ++i) f[i] *= zetan^i;

write();
write("odd terms:");
write(f);
pair[] go=fft(f,-1);
pair[] g=new pair[n];
for(int i=0; i < m; ++i) {
  g[2*i]=ge[i];
  g[2*i+1]=go[i];
}
f=fft(g,1);
// get rid of numerical error
for(int i=0; i < n; ++i) {
  if(xpart(f[i].x) < 1e-15) f[i]=(0.0,ypart(f[i]));
  if(ypart(f[i].y) < 1e-15) f[i]=(xpart(f[i]),0.0);
}
write();
write(f/n);
