include graph;

size(175,200,IgnoreAspect);

scale(Log,Log);
real[] me,e,sige,mi,i,sigi,mp,p,sigp;

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
me=a[0]; e=a[1]; sige=a[2];

file fin=input(dir+"/implicit").line();
real[][] a=fin.dimension(0,0);
a=transpose(a);
mi=a[0]; i=a[1]; sigi=a[2];

file fin=input(dir+"/pruned",check=false).line();
bool pruned=!error(fin);
if(pruned) {
  real[][] a=fin.dimension(0,0);
  a=transpose(a);
  mp=a[0]; p=a[1]; sigp=a[2];
}

colorPen[2]=heavygreen;

guide g0=scale(0.5mm)*unitcircle;
guide g1=scale(0.6mm)*polygon(3);
guide g2=scale(0.6mm)*polygon(4);

marker mark0=marker(g0,Draw(Pen(0)));
marker mark1=marker(g1,Draw(Pen(1)));
marker mark2=marker(g2,Draw(Pen(2)));

pen lp=fontsize(8pt), le=fontsize(4pt);

 // error bars:
pair[] pexp, sexp;
for(int j=0; j < me.length; ++j) {
  pexp[j]=(me[j],e[j]);
  sexp[j]=(0,e[j]-sige[j]>0?sige[j]:0);
}
errorbars(pexp,sexp,Pen(0)+le);
draw(graph(pexp),Pen(0),Label("explicit",Pen(0)+lp),mark0);
pair[] pimp, simp;
for(int j=0; j < mi.length; ++j) {
  pimp[j]=(mi[j],i[j]);
  simp[j]=(0,i[j]-sigi[j]>0?sigi[j]:0);
}
errorbars(pimp,simp,Pen(1)+le);
draw(graph(pimp),Pen(1),Label("implicit",Pen(1)+lp),mark1);
if(pruned) {
  pair[] ppru, spru;
  for(int j=0; j < mi.length; ++j) {
    ppru[j]=(mp[j],p[j]);
    spru[j]=(0,p[j]-sigp[j]>0?sigp[j]:0);
  }
  errorbars(ppru,spru,Pen(2)+le);
  draw(graph(ppru),Pen(2),Label(prunelabel,Pen(2)+lp),mark2);
}

// no error bars:
//draw(graph(me,e,e > 0),Pen(0),Label("explicit",Pen(0)+lp),mark0);
//draw(graph(mi,i,i > 0),Pen(1),Label("implicit",Pen(1)+lp),mark1);
//if(pruned) draw(graph(mp,p,p > 0),Pen(2),Label(prunelabel,Pen(2)+lp),mark2);

// fitting information; requires running rfit under R.
real[] f;
file fin=input(dir+"/implicit.p",check=false).line();
if(!error(fin)) {
  real[][] A=fin.dimension(0,0);
  real fcurve(real m) {
    real val=A[0][0]*m*log(m) +A[1][0]*m + A[2][0]*log(m) + A[3][0];
    return val;
  }

  for(int i=0; i < mi.length; ++i)
    f[i]=fcurve(mi[i]);
  // real a=min(me), b = max(me);
  // draw(graph(fcurve,a,b),Pen(1)+dashed);
  draw(graph(mi,f,f > 0),Pen(1)+dashed);
}

xaxis("$m$",BottomTop,LeftTicks);
yaxis("time (sec)",LeftRight,RightTicks);

legendlinelength=0.5cm;
legendmargin=8;
attach(legend(),point(NW),10SE);
