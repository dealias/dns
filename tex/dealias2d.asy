size(14cm,0);


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
g=box(a,2/3*b);
filldraw(shift((1/6,0))*g,gray);

z=(2,0);
g=box(a,b);
draw(shift(z)*g);
g=box(a,(1,2/3));
filldraw(shift(z)*g,gray);

g=box(a,b);
z=(4,0);
filldraw(shift(z)*g,gray);

real h=0.5;
//Arrows between boxes
draw((1,h)..(2,h),EndArrow);
label("FFTx",(1.5,h),N);
label("$\frac{2n}{3} \log \frac{2n}{3}$(?)",(1.5,h),S);
draw((3,h)..(4,h),EndArrow);
label("FFTy",(3.5,h),N);
label("$n \log \frac{2n}{3}$",(3.5,h),S);


