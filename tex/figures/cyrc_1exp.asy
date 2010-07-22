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
real dd=0.1;1.6;
real wx=0.75, wwx=1.5;
real dd=0.1;
real w=0.12;
real d=wwx+dd, mx=(wwx-wx)/2, h=dh;
pen datapen=yellow, convpen=lightgreen, labelpen=black;
path g, wg, p;
pair p1, p2, q1, q2;


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


for(int i=0; i<6; ++i){
  picture pic;
  size(pic,16cm,0);
  here=0;
  
  datapen=yellow;
  convpen=lightgreen;
  labelpen=black;
  currentpen=black;

  wg=box(a,(wwx,w));
  g=box(a,(wx,w));
  
  h=dh;
  z=(0,h);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"$\{f_n\}_{n=0}^{N-1}$",center(p),labelpen);
  p1=bottom(p);
  
  z=(d,h);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"$\{g_n\}_{n=0}^{N-1}$",center(p),labelpen);
  p2=bottom(p);
  
  iframe(i);
  h-=dh;

  z=(0,h);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"$\{f_n\}_{n=0}^{N-1}$",center(p),labelpen);
  q1=top(p);
  draw(pic,p1..q1,EndArrow);
  z+=(wx,0);
  p=shift(z)*g;
  draw(pic,p);
  label(pic,"$\{0\}_{n=0}^{N-1}$",center(p),labelpen);
  p1=bottom(p)-(wx/2,0);
  
  z=(d,h);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"$\{g_n\}_{n=0}^{N-1}$",center(p),labelpen);
  q2=top(p);
  draw(pic,p2..q2,EndArrow);
  z+=(wx,0);
  p=shift(z)*g;
  draw(pic,p);
  label(pic,"$\{0\}_{n=0}^{N-1}$",center(p),labelpen);
  p2=bottom(p)-(wx/2,0);
  

  iframe(i);
  h-=dh;

  z=(0,h);
  p=shift(z)*wg;
  filldraw(pic,p,datapen);
  label(pic,"$\{F_k\}_{k=0}^{2N-1}$",center(p),labelpen);
  q1=top(p);
  draw(pic,p1..q1,EndArrow);
  p1=bottom(p);
  
  z=(d,h);
  p=shift(z)*wg;
  filldraw(pic,p,datapen);
  label(pic,"$\{G_k\}_{k=0}^{2N-1}$",center(p),labelpen);
  q2=top(p);
  draw(pic,p2..q2,EndArrow);
  p2=bottom(p);
  
  
  iframe(i);
  h-=dh;

  z=(0,h);
  p=shift(z)*wg;
  filldraw(pic,p,convpen);
  label(pic,"$\{F_kG_k\}_{k=0}^{2N-1}$",center(p),labelpen);
  q1=top(p);
  draw(pic,p1..q1,EndArrow);
  p1=bottom(p);
  draw(pic,p2--(p2.x,0.5*(p1+q1).y)--(max(p).x,center(p).y),EndArrow);

  
  iframe(i);
  h-=dh;

  z=(0,h);
  p=shift(z)*g;
  filldraw(pic,p,convpen);
  label(pic,"$\{(f*g)_n\}_{n=0}^{N-1}$",center(p),labelpen);
  q1=top(p)+(wx/2,0);
  draw(pic,p1..q1,EndArrow);
  p1=bottom(p);

  z+=(wx,0);
  p=shift(z)*g;
  draw(pic,p);
  
  
  iframe(i);
  h-=dh;
  z=(0,h);
  p=shift(z)*g;
  filldraw(pic,p,convpen);
  label(pic,"$\{(f*g)_n\}_{n=0}^{N-1}$",center(p),labelpen);
  q1=top(p);
  draw(pic,p1..q1,EndArrow);
  
  
  shipout("cyrc_1exp"+(string) i,pic);
}
