include getparam;

import graph;
import palette;
import contour;

string[][] t={{"w","\omega"},{"vx","v_x"},{"vy","v_y"}};

string field=getstring("field","w");

// figure out how many frames there are
real[][] T;
file fin=input(run+"/t").line();
real[][] T=fin;
T=transpose(T);
int last=T[0].length-1;

int frame=getint("frame (<="+(string) last+")");
if(frame < 0) frame += last+1;
else if(frame > last) frame=last;

string name=run+"/"+field;
file fin=input(name,mode="xdr").singlereal();

real[][][] vinverted;

int nz=fin;
int ny=fin;
int nx=fin;

int pos=((3+nx*ny*nz)*frame)*4;
seek(fin,pos);
vinverted=fin.read(3);
if(eof(fin)) abort("EOF encountered on file "+name);

real[][] v;
int Ny=vinverted[0].length;

for(int j=0; j < Ny; ++j)
  v[j]=vinverted[0][Ny-1-j];

v=transpose(v);
int Nx=v[0].length;

pen[] Palette=BWRainbow2();

picture bar;
real[][] dummy;

bounds range=image(v,(0,0),(1,1),Palette,copy=false);

int Divs=10;

//real[] Cvals;
//Cvals=sequence(Divs+1)/Divs*(range.max-range.min)+range.min;
//draw(contour(v,(0,0),(1,1),Cvals,operator --));

palette(bar,Label(math(replace(field,t))),range,(0,0),
        (0.6*currentpicture.xsize,0.25cm),Bottom,
        Palette,PaletteTicks);

attach(bar.fit(),point(S),10S);
//shipout(bbox(0.25cm));
