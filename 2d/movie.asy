// Shipout each frame to separate file.

import animate;

settings.tex="pdflatex";

animation a=animation(global=false);

import graph;
import palette;
import contour;

string[][] t={{"C","C"},{"w","\omega"},{"vx","v_x"},{"vy","v_y"}};

string dir=getstring("run");
string field=getstring("field","vort");
bool global=getint("global",0) == 1;

// figure out how many frames there are
real[][] T;
file fin=input(dir+"/t").line();
real[][] T=fin;
T=transpose(T);
int count=T[0].length;

int first=getint("first frame");
int last=getint("last frame",count-1,store=false);

if(first < 0) first += count;
else if(first >= count) first=count-1;

if(last < 0) last += count;
else if(last >= count) last=count-1;

real m=infinity;
real M=-infinity;

if(global) {
  file fin=input(dir+"/"+field,mode="xdr").singlereal();
  for(int i=0; i <= last; ++i) {
    if(eof(fin))
      break;

    if(i < first) continue;

    real[][][] buf=fin.read(3);
    real[][] v=buf[0];

    m=min(m,min(v));
    M=max(M,max(v));
    buf.delete();
  }
}

file fin=input(dir+"/"+field,mode="xdr").singlereal();

pen[] Palette=BWRainbow2();

for(int i=0; i <= last; ++i) {
  picture pic;
  size(pic,20cm);

  real[][][] buf=fin.read(3);

  if(eof(fin))
    break;

  if(i < first) continue;

  real[][] v=buf[0];

  picture bar;
  bounds range=image(pic,v,global ? Range(m,M) : Full,(0,0),(1,1),Palette,
                     copy=false);
  buf.delete();

  int Divs=2;

  // real[] Cvals;
  // Cvals=sequence(Divs+1)/Divs*(range.max-range.min)+range.min;
   // draw(pic,contour(v,(0,0),(1,1),Cvals,operator --));

  palette(bar,math(replace(field,t)),range,(0,0),(0.5cm,15cm),Palette,
  	  PaletteTicks("$%+#.1f$"));
  add(pic,bar.fit(),point(pic,E),30E);
  a.add(pic);
  write("["+string(i)+"] ",none);
  purge();
}

// label(a.pdf("controls,loop"));
a.movie();
