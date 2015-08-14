include getparam;
include averages;

scale(Linear,Log);

pen p=linewidth(1);

string Etext="$E$";
string Ztext="$Z$";
string Ptext="$P$";

real[] t,E,Z,P;

real G=sqrt(2*eta)/(k0*nuH)^2;

while(nextrun()) {
  file fin=input(rundir()+"evt").line();
  real[][] a=fin.dimension(0,0);
  a=transpose(a);
  t=a[0]; E=a[1]/G^2; Z=a[2]/G^2; P=a[3];
  string runtext=" ("+run+")";
  draw(graph(t,E,E > 0),p+Pen(3*n),Etext+runtext);
  draw(graph(t,Z,Z > 0),p+Pen(3*n+1),Ztext+runtext);
}

if(n == 1) {
  currentpicture.legend[0].label=Etext;
  currentpicture.legend[1].label=Ztext;
}

xaxis("$t$",BottomTop,LeftTicks);
yaxis("Quadratic invariants",LeftRight,RightTicks(1,0));

attach(legend(),point(plain.E),20plain.E);
