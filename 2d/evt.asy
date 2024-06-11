include getparam;
include averages;

scale(Linear,Log);

pen p=linewidth(1);

string Etext="$E$";
string Ztext="$Z$";
string Ptext="$P$";

real[] t,E,Z,P,Pn;

string Pntext;
while(nextrun()) {
  Pntext="$P_"+string(nPower)+"$";
  file fin=input(rundir()+"evt").line();
  real[][] a=fin;
  a=transpose(a);
  t=a[0]; E=a[1]; Z=a[2]; P=a[3]; Pn=a[4];
  string runtext=" ("+run+")";
  draw(graph(t,E,E > 0),p+Pen(3*n),Etext+runtext);
  draw(graph(t,Z,Z > 0),p+Pen(3*n+1),Ztext+runtext);
  draw(graph(t,P,P > 0),p+Pen(3*n+2),Ptext+runtext);
  draw(graph(t,Pn,Pn > 0),p+Pen(3*n+3),Pntext+runtext);
}

if(n == 1) {
  currentpicture.legend[0].label=Etext;
  currentpicture.legend[1].label=Ztext;
  currentpicture.legend[2].label=Ptext;
  currentpicture.legend[3].label=Pntext;
}

xaxis("$t$",BottomTop,LeftTicks);
yaxis("Quadratic invariants",LeftRight,RightTicks(1,0));

attach(legend(),point(plain.E),20plain.E);
