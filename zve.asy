include getparam;
include averages;

scale(Linear,Linear);

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
  f=sqrt(2*eta/(kforce^2));
  t=a[0]; E=a[1]/f^2; Z=a[2]/f^2;
  for(int i=100; i < E.length; ++i)
    dot((E[i],Z[i]));
}

xaxis("$E$",BottomTop,LeftTicks);
yaxis("$Z$",LeftRight,RightTicks);
