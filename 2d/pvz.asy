include getparam;
include averages;

scale(Linear,Linear);
pen p=linewidth(1);

string Etext="$E$";
string Ztext="$Z$";
string Ptext="$P$";

real[] t,E,Z,P;
real f;

int k=0;
int start=getint("start");

while(nextrun()) {
  file fin=input(rundir()+"evt").line();
  real[][] a=fin;
  a=transpose(a);
  real norm=nuH^2*k0^2/eta;
  t=a[0]; E=a[1]*norm; Z=a[2]*norm; P=a[3]*norm;
  for(int i=start; i < E.length; ++i)
    dot((Z[i],P[i]),Pen(k));
  ++k;
}

xaxis("$Z$",BottomTop,LeftTicks);
yaxis("$P$",LeftRight,RightTicks);
