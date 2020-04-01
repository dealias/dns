import graph;
import palette;

size(188.80779pt,IgnoreAspect);
//size(15cm,IgnoreAspect);

int N=getint("number of timesteps:",1);
real Dt=getreal("Dt:",1);

scale(Log,Log);
pen p=linewidth(1);

real[] E=new real[N];
real[] Z=new real[N];

import utils;

real[][][] f=read(0,"F");

real f2;

real[][] f0=f[0];
real[][] f1=f[1];

pair[][] Fk=fft((pair[][]) f0,1);

int n=f0.length;
int h;

/*
for(int i=0; i < n; ++i) {
  for(int j=0; j < n; ++j) {
    if(abs(Fk[i][j]) > 1000)
      write(i,j,abs(Fk[i][j])/n^2);
  }
}
*/

for(int i=0; i < n; ++i) { 
  real[] fx=f0[i];
  real[] fy=f1[i];
  f2 += dot(fx,fx)+dot(fy,fy);
}

real nuH=0.0003;
real nuL=0.2;

write("f2=",f2);

real G=sqrt(f2)/nuH^2;

real norm=1/(G^2*nuH^2);

real[] k;

pair[][] vorticity(real[][][] U)
{
  pair[][] u=fft((pair[][]) U[0],1);
  pair[][] v=fft((pair[][]) U[1],1); 

  u=sequence(new pair[](int i) {
      pair[] ui=u[i];
      pair[] vi=v[i];
      real ki=k[i];
      return
        sequence(new pair(int j) {return I*(ki*vi[j]-k[j]*ui[j]);},n);},n);
  return fft(u,-1);
}

for(int t=0; t < N; ++t) {
  real[][][] u=read(Dt*(t+1),"U");

  real e,z;
  if(t == 0) {
    n=u[0].length;
    h=(n+1)#2;
    k=concat(sequence(0,h-1),sequence(0,h-2)-(h-1))/n^2;
  }

  pair[][] w=vorticity(u);

  real[][] u0=u[0];
  real[][] u1=u[1];
  for(int i=0; i < n; ++i) {
    real[] ux=u0[i];
    real[] uy=u1[i];
    pair[] wz=w[i];
    e += dot(ux,ux)+dot(uy,uy);
    z += dot(wz,wz).x; // Equivalent to z += sum(abs(wz)^2);
  }

  E[t]=e*norm;
  Z[t]=z*norm;
  write("["+string(t+1)+"] ",none);
}

write();
write(E,Z,sqrt(Z/E));


pen[] P=Rainbow(E.length);

for(int i=0; i < E.length; ++i) {
  frame mark;
  fill(mark,scale(0.4mm)*polygon(3),P[i]);
  add(mark,Scale((E[i],Z[i])));
}

real kfm=2.23525;
real kfp=3.60475;
  
real Emin=10^point(plain.W).x;
real Emax=10^point(plain.E).x;
real Zmax=10^point(plain.N).y;

draw(graph(new real(real E) {return E;},Emin,1),gray+p);

draw(graph(new real(real E) {return sqrt(E);},Emin,1),blue+p);

draw(graph(new real(real E) {return kfp^2*E;},
           Emin,min(10^point(plain.N).y/kfp^2,10^point(plain.E).x)),magenta);
draw(graph(new real(real E) {return kfm^2*E;},
           Emin,min(10^point(plain.N).y/kfm^2,10^point(plain.E).x)),magenta);


xaxis("$2E/(\nu G)^2$",BottomTop,LeftTicks);
yaxis("$2Z/(\nu G)^2$",LeftRight,RightTicks);
