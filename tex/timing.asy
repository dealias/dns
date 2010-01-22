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

guide g=scale(0.4mm)*unitcircle;

marker mark0=marker(g,Fill(Pen(0)));
marker mark1=marker(g,Fill(Pen(1)));
draw(graph(mp,p,p>0,Spline),Pen(0),"Padded",mark0);
draw(graph(mu,u,u>0,Spline),Pen(1),"Unpadded",mark1);

xaxis("$m$",BottomTop,LeftTicks);
yaxis("time",LeftRight,RightTicks);

attach(legend(),point(plain.E),20plain.E);

