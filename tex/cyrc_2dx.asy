size(10cm,0);


//pair center(path g) {return 0.5*(min(g)+max(g));}


void drawboxes(pair z=(0,0), pair p=(1,1))
{
  path g=box((0,0),p);
  draw(shift(z)*g);
}
path g;
pair b=(1,1), a=(0,0), z=(0,0);

g=box(a,b);
draw(shift(z)*g);
label("$f$",z+(0.5,1),N);
int n=10;
for(int i=0; i < n; ++i) {
  z=(0,i/n);
  g=box(a,(1,0.1));
  draw(shift(z)*g);
}


z=(2,0);
g=box(a,b);
draw(shift(z)*g);
label("$k_x$ even",z,SE);

g=box(a,b);
z=(3.25,0);
draw(shift(z)*g);
label("$k_x$ odd",z,SE);

real h=0.5;
//Arrows between boxes
draw((1,h)..(2,h),EndArrow);
label("FFTx",(1.5,h),N);


