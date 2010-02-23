include graph;

size(175,200,IgnoreAspect);

scale(Log,Log);
real[] mp,p,mu,u,mP,P;

string pname="cconv";
string dir;
if(pname == "conv") dir="timings1r/error.";
if(pname == "conv2") dir="timings2r/error.";
if(pname == "cconv") dir="timings1c/error.";
if(pname == "cconv2") dir="timings2c/error.";

file fin=input(dir+"explicit").line();
real[][] a=fin.dimension(0,0);
a=transpose(a);
mp=a[0]; p=a[1];

file fin=input(dir+"pruned").line();
real[][] a=fin.dimension(0,0);
a=transpose(a);
if(a.length > 1) {
  mP=a[0]; P=a[1];
}

file fin=input(dir+"implicit").line();
real[][] a=fin.dimension(0,0);
a=transpose(a);
mu=a[0]; u=a[1];

guide g0=scale(0.5mm)*unitcircle;
guide g1=scale(0.6mm)*polygon(3);
guide g2=scale(0.6mm)*polygon(4);

marker mark0=marker(g0,Draw(Pen(0)));
marker mark1=marker(g1,Draw(Pen(2)));
marker mark2=marker(g2,Draw(Pen(1)));

pen lp=fontsize(8pt);
draw(graph(mp,p,p>0),Pen(0),Label("Explicit",Pen(0)+lp),mark0);
draw(graph(mP,P,P>0),Pen(2),Label("Pruned",Pen(2)+lp),mark1);
draw(graph(mu,u,u>0),Pen(1),Label("Implicit",Pen(1)+lp),mark2);

xaxis("$m$",BottomTop,LeftTicks);
yaxis("error",LeftRight,RightTicks);

legendlinelength=0.5cm;
attach(legend(),point(NW),10SE);
