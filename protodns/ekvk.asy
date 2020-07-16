import graph;
import getrun;

size(200,150,IgnoreAspect);

int nruns=-1;

while(nextrun()) {

file in=input("zkvk").line();
real[][] a=in;
a=transpose(a);

real[] k=a[0];
real[] Zk=a[1];

scale(Log,Log);

real[] Ek=Zk/k^2;
real kstop=max(k)/sqrt(2);
draw(graph(k,Ek,k <= kstop),red);
}

xaxis("$k$",BottomTop,LeftTicks);
yaxis("$E(k)$",LeftRight,RightTicks);
