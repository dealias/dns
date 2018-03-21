include getparam;
include averages;

scale(Log,Linear);

pen p=linewidth(1);

real[] slope,kb0;
int[] index;

while(nextrun()) {
  gettime(n == 0);
  Ekavg();
  index=sequence(Ek.length-1);
  slope=log(Ek[index+1]/Ek[index])/log(k[index+1]/k[index]);
  kb0=kB[sequence(1,kB.length-2)];
  real krmax=kb0[kb0.length-1]/sqrt(2);
  draw(graph(kb0,slope,kb0 <= krmax),p+Pen(n),run);
}

ylimits(-10,0,Crop);

xaxis("$k$",BottomTop,LeftTicks);
yaxis("logarithmic slope of $E(k)$",LeftRight,RightTicks);

yequals(-3,blue+Dotted,above=false);

if(n > 1) attach(legend(),point(E),20E);
