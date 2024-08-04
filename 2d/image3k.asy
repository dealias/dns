size(0,15cm,IgnoreAspect);

import graph3;
import palette;

currentlight=light(white,specular=gray(0.7),(0,0,1));

currentprojection=orthographic(dir(30,40));
//currentprojection=orthographic(7,14,1);

//scale(Linear,Linear,Log);
usepackage("bm");

string[][] t={{"ek","$\langle\frac{|{\bm u}_{\bm k}|^2}{2}\rangle$"}};

string dir=getstring("directory","r");
string field=getstring("field","ek");

// figure out how many frames there are
real[][] T;
file fin=input(dir+"/t").line();
real[][] T=fin.dimension(0,0);
T=transpose(T);
int last=T[0].length-1;

int frame=getint("frame (<="+(string) last+")");
if (frame < 0 || frame > last) frame=last;

string name=dir+"/"+field;
file fin=input(name,mode="xdr");

int mx=fin;
int my=fin;

write(mx,my);

int N=((2mx-1)*(2my-1)-1);
int pos=(2+N)*frame*4;
seek(fin,pos);

int mx=fin;
int my=fin;

write(mx,my);

int nx=2*mx-1;
int ny=2*my-1;

real[][] v=new real[nx][ny];

for(int i=-mx+1; i < mx; ++i)
  v[mx-1+i]=new real[ny];

for(int i=-mx+1; i < mx; ++i) {
  for(int j=i <= 0 ? 1 : 0; j < my; ++j) { // start with j=1 if i <= 0
    real ek=fin;
    v[mx-1+i][my-1+j]=ek;
    v[mx-1-i][my-1-j]=ek; // Apply Hermitian symmetry
  }
}

v[mx-1][my-1]=0; // DC mode
real maxv=max(v);
if(maxv == 0) abort("Invalid spectrum.");
v[mx-1][my-1]=maxv; // override DC mode

scale(Linear,Linear,Log);

real[][] thin(real[][] v, int depth)
{
  if(depth <= 0) return v;
  int Nx=v.length;
  int Ny=v[0].length;
  int nx=Nx#2;
  int ny=Ny#2;
  real[][] V=new real[nx][];
  for(int i=0; i < nx; ++i) {
    real[] ve=v[2i];
    real[] vo=v[2i+1];
    V[i]=sequence(new real(int j) {
        return 0.25*(ve[2j]+ve[2j+1]+vo[2j]+vo[2j+1]);
      },ny);
  }
  return thin(V,depth-1);
}

real[][] V=thin(v,0);

surface s=surface(V,(-mx+1,-my+1),(mx-1,my-1));

real[] level=uniform(ScaleZ(min(V))*(1-sqrtEpsilon),
                     ScaleZ(max(V))*(1+sqrtEpsilon),256);

s.colors(palette(s.map(new real(triple v) {return find(level >= v.z);}),
                 Rainbow()));

draw(s,render(tessellate=true));

xaxis3("$k_x$",Bounds,InTicks());
yaxis3("$k_y$",Bounds,InTicks());
zaxis3(t[0][1],Bounds,InTicks);
