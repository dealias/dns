size(8cm,0);


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
label("$k_x$ even",z,SE);

g=box(a,b);
z=(1.25,0);
draw(shift(z)*g);
label("$k_x$ odd",z,SE);

frame f1;
real w=0.2;
g=box(a,(w,1));
z=(0,0);
draw(shift(z)*g,blue);
//label("$k_y$ even",z+(0,0.5),E,red);
label(f1,"$k_y$ even",z+(0,0.5));
add(rotate(90)*f1,z+(0,0.5),E);

frame f2;
z=(-1.5*w,0);
draw(shift(z)*g,blue);
label(f2,"$k_y$ odd",z+(0,0.5));
add(rotate(90)*f2,z+(0,0.5),E);


//Arrows between boxes
draw((0.05,1){NW}..{SW}(-0.15,1),EndArrow);
label("FFTy",(-0.075,1.05),N,blue);


