include getparam;
include averages;

scale(Log);

pen p=linewidth(1);

while(nextrun()) {
//  gettime(n == 0);
  gettime();
  Ekavg();
  draw(graph(k,Ek),p+Pen(2*n),texify(run));
  real A=sum(Ek);
  real B=1.0;
  real equiparition(real k) {return pi*k/(icalpha+icbeta*k^2);}
  draw(graph(equiparition,min(kc),max(kc)),p+Pen(2*n+1)+dashed,
       "$\displaystyle\frac{\pi k}{\alpha+\beta k^2}$");
}

xaxis("$k$",BottomTop,LeftTicks);
yaxis("$E(k)$",LeftRight,RightTicks);

if(currentpicture.legend.length > 1) attach(legend(),point(E),20E);
