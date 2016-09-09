include getparam;
include averages;

import palette;

scale(Linear,Linear);
pen p=linewidth(1);

string Etext="$E$";
string Ztext="$Z$";
string Ptext="$P$";

real[] t,E,Z,P;

int k=0;

while(nextrun()) {
  file fin=input(rundir()+"evt").line();
  real[][] a=fin;
  a=transpose(a);
  real norm=nuH^2*k0^2/eta;
  t=a[0]; E=a[1]*norm; Z=a[2]*norm;
  int start=getint("start",E.length#2,store=false);
  start=min(start,E.length-2);
  pen[] p=BWRainbow(E.length-start);
  frame mark;
  for(int i=start; i < E.length; ++i) {
    frame mark;
    fill(mark,scale(0.8mm)*polygon(3+k),p[i-start]);
    add(mark,(E[i],Z[i]));
  }
  draw(graph(new real(real E) {return E;},0,point(plain.E).x),blue);
  draw(graph(new real(real E) {return kforce^2*E;},
             0,min(point(plain.N).y/kforce^2,point(plain.E).x)),brown);
  draw(graph(new real(real E) {return sqrt(E);},0,(point(plain.N).y)^2),red);

  ++k;
}

xaxis("$E$",BottomTop,LeftTicks);
yaxis("$Z$",LeftRight,RightTicks);
