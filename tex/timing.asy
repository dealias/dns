include graph;

size(175,200,IgnoreAspect);

scale(Log,Log);
real[] mp,p,mu,u,mP,P;

string pname=getstring("program name");
string dir;
string prunelabel="$y$-pruned";

if(pname == "conv") dir="timings1r";
if(pname == "cconv") dir="timings1c";
if(pname == "biconv") dir="timings1b";
if(pname == "conv2") dir="timings2r";
if(pname == "cconv2") dir="timings2c";
if(pname == "biconv2") dir="timings2b";
if(pname == "cconv3") {
  dir="timings3c"; prunelabel="$xz$-pruned"; legendmargin=8;
}
  
  
file fin=input(dir+"/explicit").line();
real[][] a=fin.dimension(0,0);
a=transpose(a);
mp=a[0]; p=a[1];

file fin=input(dir+"/pruned",check=false).line();
bool pruned=!error(fin);
if(pruned) {
  real[][] a=fin.dimension(0,0);
  a=transpose(a);
  mP=a[0]; P=a[1];
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

// fitting information

file fin=input(dir+"/implicit.p").line();
real[][] A=fin.dimension(0,0);
real fcurve(real m) {
  real val=A[0][0]*m*log(m) +A[1][0]*m + A[2][0]*log(m) + A[3][0];
  return val;
}


draw(graph(mp,p,p>0),Pen(0),Label("explicit",Pen(0)+lp),mark0);

draw(graph(mu,u,u>0),Pen(1),Label("implicit",Pen(1)+lp),mark2);
real[] f;
for(int i=0; i < mu.length; ++i)
  f[i]=fcurve(mu[i]);

draw(graph(mu,f,f>0),Pen(1)+dashed);
//real a=min(mp), b = max(mp);
//draw(graph(fcurve,a,b),Pen(1)+dashed);

if(pruned) draw(graph(mP,P,P>0),Pen(2),Label(prunelabel,Pen(2)+lp),mark1);

xaxis("$m$",BottomTop,LeftTicks);
yaxis("time (sec)",LeftRight,RightTicks);

legendlinelength=0.5cm;
attach(legend(),point(NW),10SE);
