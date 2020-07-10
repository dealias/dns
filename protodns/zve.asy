size(10cm,IgnoreAspect);

import graph;
import palette;

scale(Log,Log);
pen p=linewidth(1);

real[] t,E,Z;
real G;
int k=0;
string tilde;


//file fin=input(rundir()+"evt").line();
file fin=input("ezvt").line();
real[][] a=fin;
a=transpose(a);
real fnorm=10; // Replace with sqrt(sum(Fk^2/(i*i+j*j)));
real nuH=1;
// TODO: account for nuL

G=fnorm/nuH^2;                 // Grashof number for constant forcing
tilde="";
  
write("G=",G);
real norm=G^2*nuH^2;
t=a[0]; E=2*a[1]/norm; Z=2*a[2]/norm;

int start=getint("start",a[0].length#2,store=false);
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

real E0=1e-3;
real Z0=1e-3;
draw(graph(new real(real E) {return E;},E0,1),grey);

draw(graph(new real(real E) {return sqrt(E);},Z0,1),brown);

++k;

xaxis("$2E/(\nu "+tilde+"G)^2$",BottomTop,LeftTicks);
yaxis("$2Z/(\nu "+tilde+"G)^2$",LeftRight,RightTicks);

