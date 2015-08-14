include getparam;
include averages;

scale(Linear,Log);

pen p=linewidth(1);

string Etext="$E$";
string Ztext="$Z$";
string Ptext="$P$";

real[] t,E,Z,P;
real f;

while(nextrun()) {
  file fin=input(rundir()+"evt").line();
  real[][] a=fin.dimension(0,0);
  a=transpose(a);
  f=sqrt(2*eta/kforce^2);
  t=a[0]; E=a[1]/f^2; Z=a[2]/f^2; P=a[3];
  string runtext=" ("+run+")";
  draw(graph(t,E),p+Pen(3*n),Etext+runtext);
  draw(graph(t,Z),p+Pen(3*n+1),Ztext+runtext);
}

if(n == 1) {
  currentpicture.legend[0].label=Etext;
  currentpicture.legend[1].label=Ztext;
}

xaxis("$t$",BottomTop,LeftTicks);
yaxis("Quadratic invariants",LeftRight,RightTicks(1,0));

attach(legend(),point(plain.E),20plain.E);
