size(16cm,0);


//pair center(path g) {return 0.5*(min(g)+max(g));}


void drawboxes(pair z=(0,0), pair p=(1,1))
{
  path g=box((0,0),p);
  draw(shift(z)*g);
}
path g;
pair b=(1,1), a=(0,0), z=(0,0);
real h=0.2;
real dd=0.1;
real w=0.12, wx=1.5;
real d=wx+dd;


g=box(a,(wx,w));
z=(d/2,h);
filldraw(shift(z)*g,yellow);
pair p1=z+(wx/2,0);
label("$f_n$",p1+(0,w/2));


z=(0,0);
draw(p1..(wx/2,w),EndArrow);
filldraw(shift(z)*g,yellow);
label("$F_k$",z+(wx/2,w/2));
label("$k$ even",z+(wx/2,0),S);

z=(d,0);
draw(p1..(d+wx/2,w),EndArrow);
filldraw(shift(z)*g,yellow);
label("$F_k$",z+(wx/2,w/2));
label("$k$ odd",z+(wx/2,0),S);

