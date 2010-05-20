size(18cm,0);


//pair center(path g) {return 0.5*(min(g)+max(g));}


void drawboxes(pair z=(0,0), pair p=(1,1))
{
  path g=box((0,0),p);
  draw(shift(z)*g);
}
path g;
pair b=(1,1), a=(0,0), z=(0,0);
real h=0.2;
real d=1.1;

real w=0.1;
g=box(a,(1,w));
z=(d/2,h);
draw(shift(z)*g);
label("$f_n$",z+(0.5,w/2));


z=(0,0);
draw((d/2+0.5,h)..(0.5,w),EndArrow);
//draw((0.5,w)..(0.5d+0.5,1.5w)..(d+0.5,w),EndArrow);
draw(shift(z)*g);
label("$F_k$",z+(0.5,w/2));
label("$k$ even",z+(0.5,0),S);

z=(d,0);
draw((d/2+0.5,h)..(d+0.5,w),EndArrow);
//draw((0.5,w){NE}..tension 2..{SE}(2d+0.5,w),EndArrow);
//draw((0.5,w)..(1.5d+0.5,2w)..(2d+0.5,w),EndArrow);
draw(shift(z)*g);
label("$F_k$",z+(0.5,w/2));
label("$k$ odd",z+(0.5,0),S);

