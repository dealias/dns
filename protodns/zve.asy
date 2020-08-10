size(10cm,IgnoreAspect);

import graph;
import palette;

scale(Log,Log);
pen p=linewidth(1);

real[] t,E,Z;
real G;
int k=0;
string tilde;

file fin=input("ezvt").line();

real[][] a=fin;
a=transpose(a);

int Nx=127;
int Ny=127;
int mx=(Nx+1)#2;
int my=(Ny+1)#2;

real kforce=10.0;
real deltaf=1.0;

//real fnorm;
//pair F=(0.3,0);
//real F2=abs(F)^2;
//for(int i=-mx+1; i < mx; ++i) {
//  int i2=i*i;
//  real halfdeltaf=0.5*deltaf;
//  for(int j=(i <= 0 ? 1 : 0); j < my; ++j) {
//    real k=sqrt(i2+j*j);
//    fnorm += abs(k-kforce) < halfdeltaf ? F2/(k*k) : 0.0;
//  }
//}

// Account for reality condition
//fnorm=sqrt(2*fnorm);

real nuH=0.003;
real nuL=0.2;

real eta=1.0;
real eps=1.0/(kforce)^2;
G=sqrt(eps*(nuH+nuL))/nuH^2;                // Grashof number for stochastic forcing
tilde="";

write("G=",G);
real norm=G^2*nuH^2;
t=a[0]; E=2*a[1]/norm; Z=2*a[2]/norm;

int start=getint("start",0,store=false);
int end=getint("end",a[0].length,store=false);
t=t[start:end];
real t0=t[0];
real tmax=t[t.length-1];
write("t:",t0,tmax);
E=E[start:end];
Z=Z[start:end];
real incr=(E.length-1)/tmax;
pen[] p=Rainbow(E.length);

for(int i=0; i < E.length; ++i) {
  frame mark;
  fill(mark,scale(0.4mm)*polygon(3+k),p[round(t[i]*incr)]);
  add(mark,Scale((E[i],Z[i])));
}

real Emin=10^point(plain.W).x;
real Emax=10^point(plain.E).x;
real Zmin=10^point(plain.S).y;
real Zmax=10^point(plain.N).y;

draw(graph(new real(real E) {return E;},Emin,1),grey);

draw(graph(new real(real E) {return sqrt(E);},Zmin,1),brown);

++k;

xaxis("$2E/(\nu "+tilde+"G)^2$",BottomTop,LeftTicks);
yaxis("$2Z/(\nu "+tilde+"G)^2$",LeftRight,RightTicks);