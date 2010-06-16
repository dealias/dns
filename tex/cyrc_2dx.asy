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


  z=(0,dd);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"$f$",center(p),labelpen);
  pair f1=(center(p).x,min(p).y);
  
  real ddh=0.1dh+2dd;
  
  z=(ddh,dd);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"$g$",center(p),labelpen);
  pair g1=(center(p).x,min(p).y);

  iframe(i);
  
  z=(0,0);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"${\rm FFT}_x\{f\}$",center(p),N);
  label(pic,"$k_x$ \rm{even}",center(p),S);
  pair ff1=(center(p).x,max(p).y);
  draw(pic,f1..ff1,EndArrow);
  
  z=(dh,0);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"${\rm FFT}_x\{f\}$",center(p),N);
  label(pic,"$k_x$ \rm{odd}",center(p),S);
  pair ff2=(center(p).x,max(p).y);
  draw(pic,f1..ff2,EndArrow);
  
  z=(ddh,0);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"${\rm FFT}_x\{g\}$",center(p),N);
  label(pic,"$k_x$ \rm{even}",center(p),S);
  pair ff1=(center(p).x,max(p).y);
  draw(pic,g1..ff1,EndArrow);
  
  z=(ddh+dh,0);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"${\rm FFT}_x\{g\}$",center(p),N);
  label(pic,"$k_x$ \rm{odd}",center(p),S);
  pair ff2=(center(p).x,max(p).y);
  draw(pic,g1..ff2,EndArrow);

  shipout("cyrc_2dx"+(string) i,pic);
}
