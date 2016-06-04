include getparam;
include averages;

scale(Log);

pen p=linewidth(1);

while(nextrun()) {
  //  gettime(n == 0);
  gettime();
  Ekavg();

  real krmax=spectrum != 3 ? kb[kb.length-1]/sqrt(2) : kb[kb.length-1]+1;
 
  draw(graph(k,Ek,k <= krmax),p+Pen(n),texify(run));
  //draw(graph(k,Ek,k <= krmax),p+Pen(n),texify(run),marker(scale(0.8mm)*polygon(3+n)));
}

xaxis("$k$",BottomTop,LeftTicks);
yaxis("$E(k)$",LeftRight,RightTicks);

if(currentpicture.legend.length > 1) attach(legend(),point(E),20E);
