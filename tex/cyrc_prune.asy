size(10cm,0);


pair center(path g) {return 0.5*(min(g)+max(g));}


void drawboxes(pair z=(0,0), pair p=(1,1))
{
  path g=box((0,0),p);
  draw(shift(z)*g);
}
path g, p;
pair b=(1,1), a=(0,0), z=(0,0);

g=box(a,b);
draw(shift(z)*g);
g=box(a,0.5*b);
filldraw(shift(z)*g,yellow);
g=box(a,(1,0.5));
p=shift(z)*g;
draw(p);
pair p1=(max(p).x,center(p).y);


z=(2,0);
g=box(a,(1,0.5));
p=shift(z)*g;
filldraw(shift(z)*g,yellow);
pair p2=(min(p).x,center(p).y);
draw(p1..p2,EndArrow);
label("x-FFT",(p1+p2)/2,N);
g=box(a,b);
p=shift(z)*g;
draw(p);
p1=(max(p).x,center(p).y);


z=(4,0);
g=box(a,b);
p=shift(z)*g;
filldraw(p,yellow);
pair p2=(min(p).x,center(p).y);
draw(p1..p2,EndArrow);

label("y-FFT",(p1+p2)/2,N);

/*
int n=10;
for(int i=0; i < n-1; ++i) {
  z=(2+i/n,0);
  g=box(a,(0.1,1));
  draw(shift(z)*g,blue);
}
*/



/*
g=box(a,b);
draw(shift(z)*g);
int n=10;
for(int i=0; i < n-1; ++i) {
  z=(0,0.5*i/n);
  g=box(a,(1,0.1));
  draw(shift(z)*g,red);
}
*/


