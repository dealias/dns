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
  real[][] a=fin.dimension(0,0);
  a=transpose(a);
  f=sqrt(2*eta/(kforce^2));
  t=a[0]; E=a[1]/f^2; Z=a[2]/f^2;
  for(int i=start; i < E.length; ++i)
    dot((E[i],Z[i]),Pen(k));
  ++k;
}

xaxis("$E$",BottomTop,LeftTicks);
yaxis("$Z$",LeftRight,RightTicks);
