include graph;

size(175,200,IgnoreAspect);

scale(Log,Log);
real[] mp,p,mu,u,mP,P;

string pname=getstring("program name");
string dir;
bool prune=false;
if(pname == "conv") dir="timings1r";
if(pname == "cconv") dir="timings1c";
if(pname == "conv2") {dir="timings2r";prune=true;}
if(pname == "cconv2") {dir="timings2c";prune=true;}
if(pname == "biconv") dir="timings1b";
if(pname == "biconv2") dir="timings2b";
  
file fin=input(dir+"/explicit").line();
real[][] a=fin.dimension(0,0);
a=transpose(a);
mp=a[0]; p=a[1];

if(prune) {
  file fin=input(dir+"/pruned").line();
  real[][] a=fin.dimension(0,0);
  a=transpose(a);
  if(a.length > 1) {
    mP=a[0]; P=a[1];
  }
}

file fin=input(dir+"/implicit").line();
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
if(prune) draw(graph(mP,P,P>0),Pen(2),Label("Pruned",Pen(2)+lp),mark1);
draw(graph(mu,u,u>0),Pen(1),Label("Implicit",Pen(1)+lp),mark2);

xaxis("$m$",BottomTop,LeftTicks);
yaxis("time (sec)",LeftRight,RightTicks);

legendlinelength=0.5cm;
attach(legend(),point(NW),10SE);
