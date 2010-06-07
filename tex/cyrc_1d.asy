size(16cm,0);


pair center(path g) {return 0.5*(min(g)+max(g));}


void drawboxes(pair z=(0,0), pair p=(1,1))
{
  path g=box((0,0),p);
  draw(shift(z)*g);
}
path g, p;
pair b=(1,1), a=(0,0), z=(0,0);
real h=0.3;
real dd=0.1;
real w=0.12, wx=1.5;
real d=wx+dd;


g=box(a,(wx,w));
z=(d/2,h);
p=shift(z)*g;
filldraw(p,yellow);
pair p1=center(p);
label("$f_n$",center(p));

z=(0,0);
p=shift(z)*g;
filldraw(p,yellow);
draw(p1-(0,w/2)..center(p)+(0,w/2),EndArrow);
label("$F_k$, $k$ even",center(p));

z=(d,0);
p=shift(z)*g;
filldraw(p,yellow);
draw(p1-(0,w/2)..center(p)+(0,w/2),EndArrow);
label("$F_k$, $k$ odd",center(p));

