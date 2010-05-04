include getparam;
include averages;

scale(Log);

pen p=linewidth(1);



while(nextrun()) {
//  gettime(n == 0);
  gettime();
  Ekavg();
  draw(graph(k,Ek),p+Pen(n),texify(run));
  real A=sum(Ek);
  real B=1.0;
  //  draw(graph(k,k^2/(A+B*k^2)),p+Pen(n)+dashed);
  real equiparition(real k) {  return pi*k/(icalpha+icbeta*k*k);}
  draw(graph(equiparition,min(kc),max(kc)),p+Pen(n)+dashed);
}

xaxis("$k$",BottomTop,LeftTicks);
yaxis("$E(k)$",LeftRight,RightTicks);

if(currentpicture.legend.length > 1) attach(legend(),point(E),20E);
