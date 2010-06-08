size(16cm,0);


//pair center(path g) {return 0.5*(min(g)+max(g));}


void drawboxes(pair z=(0,0), pair p=(1,1))
{
  path g=box((0,0),p);
  draw(shift(z)*g);
}
path g, wg;
pair b=(1,1), a=(0,0), z=(0,0);
real dh=0.25;
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


for(int i=0; i<3; ++i){
  picture pic;
  size(pic,8cm,0);
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
  label(pic,"$\{f_n\}_{n=0}^{N-1}$",p1+(0,w/2),labelpen);


  
  h-=dh;

  iframe(i);

  z=(0,h);
  filldraw(pic,shift(z)*g,datapen);
  pair p1e=z+(wx/2,0)+(0,w);
  label(pic,"$\{F_k\}_{n=0}^{N-1}, k\, \rm{even}$",p1e-(0,w/2),labelpen);
  draw(pic,p1..p1e,EndArrow);

  iframe(i);
  
  z=(wx+dd,h);
  filldraw(pic,shift(z)*g,datapen);
  pair p1o=z+(wx/2,0)+(0,w);
  label(pic,"$\{F_k\}_{n=0}^{N-1}, k\, \rm{odd}$",p1o-(0,w/2),labelpen);
  draw(pic,p1..p1o,EndArrow);
  
  
  
  shipout("cyrc_1d"+(string) i,pic);
}
