size(0,15cm,IgnoreAspect);

import graph3;
import palette;

currentlight=Viewport;

currentprojection=orthographic(dir(72,40));

usepackage("bm");

string dir=getstring("directory","");
string field=getstring("field","angle");

// figure out how many frames there are
real[][] T;
file fin=input(dir+"/t").line();
real[][] T=fin.dimension(0,0);
T=transpose(T);
int last=T[0].length-1;

int frame=getint("frame (<= "+(string) last+")");
if (frame < 0 || frame > last) frame=last;

string name=dir+"/"+field;
file fin=input(name,mode="xdr").singlereal();

int nx=fin;
int ny=fin;

int pos=(2+nx*ny)*frame*4;
seek(fin,pos);

int ny=fin;
int nx=fin;
write(nx,ny);

real[][] v=fin.dimension(nx,ny);

if(eof(fin)) abort("EOF encountered on file "+name);

int Nx=1024;
int Ny=Nx;

real[][] V=v[0:Nx];

for(int j=0; j < V.length; ++j)
  V[j]=V[j][0:Ny];

write(V.length);
write(V[0].length);

real maxV=max(V);
real minV=min(V);

int Ncolors=256;

real[] level=uniform(ScaleZ(minV)*(1-sqrtEpsilon),
                     ScaleZ(maxV)*(1+sqrtEpsilon),Ncolors);

surface s=surface(V,(0,0),(Nx,Ny),Spline);
s.colors(palette(s.map(new real(triple v) {return find(level >= v.z);}),
                 BWRainbow2()));
draw(s,render(tessellate=true));

usepackage("bm");
texpreamble("
\def\v{\bm}
\def\vu{{\v u}}
\def\B{{\cal B}}
");

Label zlabel;
if(field == "w")
  zlabel="$\omega$";
if(field == "angle")
  zlabel="$\langle\theta\rangle$";
if(field == "triplet")
  zlabel=rotate(90)*"$\frac{(\B(\vu,\vu),A^2\vu)}{2\nu P_2}$";

xaxis3("$x$",Bounds,InTicks(endlabel=false));
yaxis3("$y$",Bounds,InTicks);
zaxis3(zlabel,Bounds,InTicks);
