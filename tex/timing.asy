include graph;
size(200,150,IgnoreAspect);

scale(Log,Log);
real[] mp,p,mu,u;

file fin=input("padded").line();
real[][] a=fin.dimension(0,0);
a=transpose(a);
mp=a[0]; p=a[1];

file fin=input("unpadded").line();
real[][] a=fin.dimension(0,0);
a=transpose(a);
mu=a[0]; u=a[1];

marker mark=marker(scale(0.8mm)*unitcircle);
draw(graph(mp,p,p>0),Pen(0),"Padded",mark);
draw(graph(mu,u,u>0),Pen(1),"Unpadded",mark);

xaxis("$m$",BottomTop,LeftTicks);
yaxis("time",LeftRight,RightTicks);

attach(legend(),point(plain.E),20plain.E);

