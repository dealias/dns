size(10cm,0);


//pair center(path g) {return 0.5*(min(g)+max(g));}


void drawboxes(pair z=(0,0), pair p=(1,1))
{
  path g=box((0,0),p);
  draw(shift(z)*g);
}
path g;
pair b=(1,1), a=(0,0), z=(0,0);




real l=0.1;
real h=0;
real h2=-0.4;
real x2=1.1;
g=box(a,(1,l));

z=(0,h);
draw(shift(z)*g,red);
label("$F_k$: $k_y$ even",z+(0.5,0.1),N);

z=(x2,h);
draw(shift(z)*g,red);
label("$F_k$: $k_y$ odd",z+(0.5,0.1),N);

label("$\times$",(1+(x2-1)/2,l/2+h2/2));

z=(0,h2);
draw(shift(z)*g,red);
label("$G_k$: $k_y$ even",z+(0.5,0.1),N);

z=(x2,h2);
draw(shift(z)*g,red);
label("$G_k$: $k_y$ odd",z+(0.5,0.1),N);

label("$=$",(1+(x2-1)/2,l/2+h2+h2/2));

z=(0,2*h2);
draw(shift(z)*g,blue);
label("$F_k G_k$: $k_y$ even",z+(0.5,0.1),N);
z=(x2,2*h2);
draw(shift(z)*g,blue);
label("$F_k G_k$: $k_y$ odd",z+(0.5,0.1),N);
