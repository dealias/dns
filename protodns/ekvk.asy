import graph;

size(200,150,IgnoreAspect);

file in=input("Zk.dat").line();
real[][] a=in.dimension(0,0);
a=transpose(a);

real[] k=a[0];
real[] Zk=a[1];

scale(Log,Log);

real[] Ek=Zk/k^2;
draw(graph(k,Ek),red);

xaxis("$k$",BottomTop,LeftTicks);
yaxis("$E(k)$",LeftRight,RightTicks);
