include getparam;
include averages;

scale(Linear,Linear);
pen p=linewidth(1);

string Etext="$E$";
string Ztext="$Z$";
string Ptext="$P$";

real[] t,E,Z,P;

int k=0;
int start=getint("start");

while(nextrun()) {
  file fin=input(rundir()+"evt").line();
  real[][] a=fin;
  a=transpose(a);
  real norm=nuH^2*k0^2/eta;
  t=a[0]; E=a[1]*norm; Z=a[2]*norm;
  dot((E[start],Z[start]),Pen(k));
  for(int i=start+1; i < E.length; ++i) {
    dot((E[i],Z[i]),Pen(k));
    //          draw((E[i-1],Z[i-1])--(E[i],Z[i]),blue,Arrow);
  }
  draw(graph(new real(real E) {return E;},0,point(plain.E).x),blue);
  draw(graph(new real(real E) {return sqrt(E);},0,(point(plain.N).y)^2),blue);

  ++k;
}

xaxis("$E$",BottomTop,LeftTicks);
yaxis("$Z$",LeftRight,RightTicks);
