// Shipout each frame to separate file.

import animate;

bool transpose=false;

settings.tex="pdflatex";

animation a=animation(global=false);

import graph;
import palette;
import contour;

string[][] t={{"C","C"},{"w","\omega"},{"vx","v_x"},{"vy","v_y"}};

string dir=getstring("directory","r");
string field=getstring("field","vort");

// figure out how many frames there are
real[][] T;
file fin=input(dir+"/t").line();
real[][] T=fin.dimension(0,0);
T=transpose(T);
int count=T[0].length;

int first=getint("first frame");
int last=getint("last frame");

if(first < 0) first += count;
else if(first >= count) first=count-1;

if(last < 0) last += count;
else if(last >= count) last=count-1;

file fin=input(dir+"/"+field,check=true,mode="xdr").singlereal();

bounds range;
pen[] Palette=BWRainbow2();

for(int i=0; i <= last; ++i) {
  picture pic;
  real[][][] buf;
  size(pic,300);
 
  if(transpose) {
    int nx=fin;
    int ny=fin;
    int nz=fin;
    buf=fin.dimension(nz,ny,nx);
  } else
    buf=fin.read(3);

  if(eof(fin))
    break;

  if(i < first) continue;

  real[][] v=buf[0];

  picture bar;
  bounds thisrange=image(pic,v,(0,0),(1,1),Palette,transpose=false,copy=false);

  if(i == first) range=thisrange;

  int Divs=2;

  // real[] Cvals;
  // Cvals=sequence(Divs+1)/Divs*(range.max-range.min)+range.min;
  // draw(pic,contour(v,(0,0),(1,1),Cvals,operator --));

  palette(bar,math(replace(field,t)),range,(0,0),(0.5cm,8cm),Palette,
  	  PaletteTicks("$%+#.1f$"));
  add(pic,bar.fit(),point(pic,E),30E);
  a.add(pic);
  write("["+string(i)+"] ",none);
  purge();
}

// label(a.pdf("controls,loop"));
a.movie();
