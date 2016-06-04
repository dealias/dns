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
  f=sqrt(2*eta/kforce^2);
  t=a[0]; E=a[1]/f^2; Z=a[2]/f^2;
  dot((E[start],Z[start]),Pen(k));
  for(int i=start+1; i < E.length; ++i) {
    dot((E[i],Z[i]),Pen(k));
    draw((E[i-1],Z[i-1])--(E[i],Z[i]),blue,Arrow);
  }

  ++k;
}

xaxis("$E$",BottomTop,LeftTicks);
yaxis("$Z$",LeftRight,RightTicks);
