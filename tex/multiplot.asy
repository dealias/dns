include graph;
size(200,150,IgnoreAspect);

scale(Log,Log);

string fdata=getstring("filename",default="");
int n=0;
guide g=scale(0.4mm)*unitcircle;
while(fdata != "") {
  file fin=input(fdata).line();
  real[][] a=fin.dimension(0,0);
  a=transpose(a);
  real[] a0=a[0];
  real[] a1=a[1];
  marker mark0=marker(g,Fill(Pen(n)));
  draw(graph(a0,a1,a1>0),Pen(n),fdata,mark0);
  ++n;
  fdata=getstring("next filename");
  if(fdata == "no") fdata="";
}
xaxis("$m$",BottomTop,LeftTicks);
yaxis("time",LeftRight,RightTicks);

attach(legend(),point(plain.E),20plain.E);
