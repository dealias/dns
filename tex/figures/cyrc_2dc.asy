size(11cm,0);


pair center(path g) {return 0.5*(min(g)+max(g));}


void drawboxes(pair z=(0,0), pair p=(1,1))
{
  path g=box((0,0),p);
  draw(shift(z)*g);
}
path g, p;
pair b=(1,1), a=(0,0), z=(0,0);
real wx=0.01;
real dd=1.2*wx;
real dh=1.1*wx;
real ddh=0.1dh+2dd;
g=box(a,(wx,wx));
//draw(shift(z)*g);

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


for(int i=0; i<2; ++i){
  picture pic;
  size(pic,11cm,0);
  here=0;
  
  datapen=yellow;
  convpen=lightgreen;
  labelpen=black;
  currentpen=black;
  
  z=(0,0);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"${\rm FFT}_x\{f\}$",center(p),N,labelpen);
  label(pic,"$k_x$ \rm{even}",center(p),S,labelpen);
  pair ff1=(center(p).x,min(p).y);
  
  z=(dh,0);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"${\rm FFT}_x\{f\}$",center(p),N,labelpen);
  label(pic,"$k_x$ \rm{odd}",center(p),S,labelpen);
  pair ff2=(center(p).x,min(p).y);
  
  z=(ddh,0);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"${\rm FFT}_x\{g\}$",center(p),N,labelpen);
  label(pic,"$k_x$ \rm{even}",center(p),S,labelpen);
  pair gg1=(center(p).x,min(p).y);
  
  
  z=(ddh+dh,0);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"${\rm FFT}_x\{g\}$",center(p),N,labelpen);
  label(pic,"$k_x$ \rm{odd}",center(p),S,labelpen);
  pair gg2=(center(p).x,min(p).y);

  iframe(i);
  
  z=(0,-dd);
  p=shift(z)*g;
  filldraw(pic,p,convpen);
  label(pic,"${\rm FFT}_x\{f*g\}$",center(p),N,labelpen);
  label(pic,"$k_x$ \rm{even}",center(p),S,labelpen);
  pair f1=(center(p).x,max(p).y);
  draw(pic,ff1..f1,EndArrow);
  draw(pic,gg1..f1,EndArrow);
  
  z=(dh,-dd);
  p=shift(z)*g;
  filldraw(pic,p,convpen);
  label(pic,"${\rm FFT}_x\{f*g\}$",center(p),N,labelpen);
  label(pic,"$k_x$ \rm{odd}",center(p),S,labelpen);
  pair g1=(center(p).x,max(p).y);
  draw(pic,ff2..g1,EndArrow);
  draw(pic,gg2..g1,EndArrow);


  shipout("cyrc_2dc"+(string) i,pic);
}
