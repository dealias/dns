size(0,4cm);


pair center(path g) {return 0.5*(min(g)+max(g));}
pair bottom(path g) {return (center(g).x,min(g).y);}


void drawboxes(pair z=(0,0), pair p=(1,1))
{
  path g=box((0,0),p);
  draw(shift(z)*g);
}
path g, p;
pair b=(1,1/2), a=(0,0), z=(0,0);
real h=0.5;

g=box(a,b);
draw(shift(z)*g);

g=box(a,2/3*b);
z=(1/6,0);
p=shift(z)*g;
filldraw(p,yellow);

dot(bottom(p),red);
