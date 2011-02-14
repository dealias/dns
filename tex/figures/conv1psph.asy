size(16cm,0);


pair center(path g) {return 0.5*(min(g)+max(g));}
pair bottom(path g) {return (center(g).x,min(g).y);}
pair top(path g) {return (center(g).x,max(g).y);}


void drawboxes(pair z=(0,0), pair p=(1,1))
{
  path g=box((0,0),p);
  draw(shift(z)*g);
}

pair b=(1,1), a=(0,0), z=(0,0);
real dh=0.25;
real dd=0.1;
real wx=0.75, wwx=1.5;
real dd=0.1;
real w=0.12;
real d=wx+dd, mx=(wwx-wx)/2, h=dh;
pen datapen=yellow, convpen=lightgreen, labelpen=black;
path g, p;
pair p1, p2, p3, p4, q1, q2, q2, q3, q4;


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


for(int i=0; i<5; ++i){
  picture pic;
  size(pic,16cm,0);
  here=0;
  
  datapen=yellow;
  convpen=lightgreen;
  labelpen=black;
  currentpen=black;


  g=box(a,(wx,w));
  
  h=dh;

  z=(0,h);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"$\{F_k\}_{k=0}^{N-1}$",center(p),labelpen);
  p1=bottom(p);
  p3=bottom(p);
  
  z=(d,h);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"$\{G_k\}_{k=0}^{N-1}$",center(p),labelpen);
  p2=bottom(p);
  p4=bottom(p);


  iframe(i);
  h-=dh;

  z=(0,h);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"$\{f_n\}$",center(p),labelpen);
  q1=top(p);
  draw(pic,p1..q1,EndArrow);
  p1=bottom(p);
  
  z+=(d,0);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"$\{g_n\}$",center(p),labelpen);
  q2=top(p);
  draw(pic,p2..q2,EndArrow);
  p2=bottom(p);
  
  z+=(d,0);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"$\{f^\Delta_k\}$",center(p),labelpen);
  q3=top(p);
  draw(pic,p3..q3,EndArrow);
  p3=bottom(p);
  
  z+=(d,0);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"$\{g^\Delta_k\}$",center(p),labelpen);
  q4=top(p);
  draw(pic,p4..q4,EndArrow);
  p4=bottom(p);


  iframe(i);
  h-=dh;

  z=(0,h);
  p=shift(z)*g;
  filldraw(pic,p,convpen);
  label(pic,"$\{f_n g_n\}$",center(p),labelpen);
  q1=top(p);
  draw(pic,p1..q1,EndArrow);
  draw(pic,p2..q1,EndArrow);
  p1=bottom(p);
    
  z+=(2d,0);
  p=shift(z)*g;
  filldraw(pic,p,convpen);
  label(pic,"$\{f^\Delta_k f^\Delta_k\}$",center(p),labelpen);
  q3=top(p);
  draw(pic,p3..q3,EndArrow);
  draw(pic,p4..q3,EndArrow);
  p3=bottom(p);

  iframe(i);
  h-=dh;

  z=(0,h);
  p=shift(z)*g;
  filldraw(pic,p,convpen);
  label(pic,"$\{F *_N G\}$",center(p),labelpen);
  q1=top(p);
  draw(pic,p1..q1,EndArrow);
  p1=bottom(p);
    
  z+=(2d,0);
  p=shift(z)*g;
  filldraw(pic,p,convpen);
  label(pic,"$\{F *_\Delta G\}$",center(p),labelpen);
  q3=top(p);
  draw(pic,p3..q3,EndArrow);
  p3=bottom(p);


  iframe(i);
  h-=dh;

  z=(0,h);
  p=shift(z)*g;
  filldraw(pic,p,convpen);
  label(pic,"$\{F*G\}$",center(p),labelpen);
  q1=top(p);
  draw(pic,p1..q1,EndArrow);
  draw(pic,p3..q1,EndArrow);
  p1=bottom(p);

  
  shipout("conv1psph"+(string) i,pic);
}
