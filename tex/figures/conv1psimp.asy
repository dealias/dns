size(16cm,0);


//pair center(path g) {return 0.5*(min(g)+max(g));}


void drawboxes(pair z=(0,0), pair p=(1,1))
{
  path g=box((0,0),p);
  draw(shift(z)*g);
}
path g, wg;
pair b=(1,1), a=(0,0), z=(0,0);
real dh=0.35;
real dd=0.1;1.6;
real wx=0.75, wwx=1.5;
real dd=0.1;
real w=0.12;
real d=wwx+dd, mx=(wwx-wx)/2, h=dh;
pen datapen=yellow, convpen=lightgreen, labelpen=black;


int here=0;
void iframe(int i) {
  if (here == i) {
    datapen=invisible;
    convpen=invisible;
    labelpen=invisible;
    currentpen=invisible;
  }
  ++here;
}


for(int i=0; i<4; ++i){
  picture pic;
  size(pic,16cm,0);
  here=0;
  
  datapen=yellow;
  convpen=lightgreen;
  labelpen=black;
  currentpen=black;

  h=dh;
  g=box(a,(wx,w));
  //  z=((wx+dd)/2,h);
  z=(0,h);
  filldraw(pic,shift(z)*g,datapen);
  pair p1=z+(wx/2,0);
  label(pic,"$\{F_k\}_{k=0}^{N-1}$",p1+(0,w/2),labelpen);

  //  z=(d+(wx+dd)/2,h);
  z=(d+dd,h);
  filldraw(pic,shift(z)*g,datapen);
  pair p2=z+(wx/2,0);
  label(pic,"$\{G_k\}_{k=0}^{N-1}$",p2+(0,w/2),labelpen);

  wg=box(a,(wwx,w));

  iframe(i);
  h-=dh;
  
  z=(0,h);
  filldraw(pic,shift(z)*g,datapen);
  pair p1e=z+(wx/2,0)+(0,w);
  label(pic,"$\{f_n\}_{n=0}^{N-1}, n\, \rm{even}$",p1e-(0,w/2),labelpen);
  draw(pic,p1..p1e,EndArrow);
  
  z=(wx+dd,h);
  filldraw(pic,shift(z)*g,datapen);
  pair p1o=z+(wx/2,0)+(0,w);
  label(pic,"$\{f_n\}_{n=0}^{N-1}, n\, \rm{odd}$",p1o-(0,w/2),labelpen);
  draw(pic,p1..p1o,EndArrow);
  
  z=(2wx+2dd,h);
  filldraw(pic,shift(z)*g,datapen);
  pair p2e=z+(wx/2,0)+(0,w);
  label(pic,"$\{g_n\}_{n=0}^{N-1}, n\, \rm{even}$",p2e-(0,w/2),labelpen);
  draw(pic,p2..p2e,EndArrow);
  
  z=(3wx+3dd,h);
  filldraw(pic,shift(z)*g,datapen);
  pair p2o=z+(wx/2,0)+(0,w);
  label(pic,"$\{g_n\}_{n=0}^{N-1}, n\, \rm{odd}$",p2o-(0,w/2),labelpen);
  draw(pic,p2..p2o,EndArrow);

  
  iframe(i);
  h-=dh;

  //  z=(wx+dd,h);
  z=(0,h);
  filldraw(pic,shift(z)*g,convpen);
  pair Pe=z+(wx/2,0)+(0,w);
  label(pic,"$\{f_n g_n\}_{n=0}^{N-1}, n\, \rm{even}$",Pe-(0,w/2),labelpen);
  draw(pic,p2e-(0,w)..Pe,EndArrow);
  draw(pic,p1e-(0,w)..Pe,EndArrow);
  
  //z=(2wx+2dd,h);
  z=(wx+dd,h);
  filldraw(pic,shift(z)*g,convpen);
  pair Po=z+(wx/2,0)+(0,w);
  label(pic,"$\{f_n g_n\}_{k=0}^{N-1}, n\, \rm{odd}$",Po-(0,w/2),labelpen);
  draw(pic,p2o-(0,w)..Po,EndArrow);
  draw(pic,p1o-(0,w)..Po,EndArrow);
    
  iframe(i);
  h-=dh;
  //  z=(1.5wx+1.5dd,h);
  z=(0,h);
  filldraw(pic,shift(z)*g,convpen);
  pair plast=z+(wx/2,0)+(0,w);
  label(pic,"$\{(F*G)_k\}_{k=0}^{N-1}$",plast-(0,w/2),labelpen);
  draw(pic,Pe-(0,w)..plast,EndArrow);
  draw(pic,Po-(0,w)..plast,EndArrow);
  
  
  shipout("conv1psimp"+(string) i,pic);
}
