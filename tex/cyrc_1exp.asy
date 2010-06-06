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


for(int i=0; i<6; ++i){
  picture pic;
  size(pic,16cm,0);
  here=0;
  
  datapen=yellow;
  convpen=lightgreen;
  labelpen=black;
  currentpen=black;

  h=dh;
  g=box(a,(wx,w));
  z=(mx,h);
  filldraw(pic,shift(z)*g,datapen);
  label(pic,"$\{f_n\}_{n=0}^{N-1}$",z+(wx/2,w/2),labelpen);

  
  z=(mx+d,h);
  filldraw(pic,shift(z)*g,datapen);
  label(pic,"$\{g_n\}_{n=0}^{N-1}$",z+(wx/2,w/2),labelpen);


  iframe(i);
  h-=dh;
  wg=box(a,(wwx,w));
  z=(0,h);
  
  draw(pic,shift(z)*wg);
  filldraw(pic,shift(z)*g,datapen);
  label(pic,"$\{f_n\}_{n=0}^{N-1}$",z+(wx/2,w/2),labelpen);
  label(pic,"$\{0\}_{n=0}^{N-1}$",z+(wx+mx,w/2),labelpen);
  draw(pic,z+(wwx/2,dh)..z+(wwx/2,0)+(0,w),EndArrow);
  // draw(pic,z+(wwx/2,0)..z+(wwx/2,0)+(0,w-dh),EndArrow);

  z=(d,h);
  draw(pic,shift(z)*wg);
  filldraw(pic,shift(z)*g,datapen);
  label(pic,"$\{g_n\}_{n=0}^{N-1}$",z+(wx/2,w/2),labelpen);
  label(pic,"$\{0\}_{n=0}^{N-1}$",z+(wx+mx,w/2),labelpen);
  draw(pic,z+(wwx/2,dh)..z+(wwx/2,0)+(0,w),EndArrow);
  //  draw(pic,z+(wwx/2,0)..z+(wwx/2,0)+(0,w-dh),EndArrow);

  iframe(i);
  h-=dh;
  z=(0,h);
  filldraw(pic,shift(z)*wg,datapen);
  label(pic,"$\{F_k\}_{k=0}^{N-1}$",z+(wwx/2,w/2),labelpen);
  draw(pic,z+(wwx/2,dh)..z+(wwx/2,0)+(0,w),EndArrow);

  
  z=(d,h);
  filldraw(pic,shift(z)*wg,datapen);
  label(pic,"$\{G_k\}_{k=0}^{2N-1}$",z+(wwx/2,w/2),labelpen);
  draw(pic,z+(wwx/2,dh)..z+(wwx/2,0)+(0,w),EndArrow);
  
  
  iframe(i);
  h-=dh;
  z=((wwx+dd)/2,h);

  pair dest=(wwx+dd/2,w+h);
  draw(pic,(wwx/2,h+dh)..dest,EndArrow);
  draw(pic,(1.5*wwx+dd,h+dh)..dest,EndArrow);


  //z=(0,h);
  filldraw(pic,shift(z)*wg,convpen);
  label(pic,"$\{F_kG_k\}_{k=0}^{2N-1}$",z+(wwx/2,w/2),labelpen);

  iframe(i);
  draw(pic,z+(wwx/2,0)..z+(wwx/2,0)+(0,w-dh),EndArrow);
  h-=dh;
  z=((wwx+dd)/2,h);
  //z=(0,h);
  draw(pic,shift(z)*wg);
  filldraw(pic,shift(z)*g,convpen);
  label(pic,"$\{f*g_n\}_{n=0}^{N-1}$",z+(wx/2,w/2),labelpen);
  label(pic,"$\{0\}_{n=0}^{N-1}$",z+(wx+mx,w/2),labelpen);


  iframe(i);
  draw(pic,z+(wx/2,0)..z+(wwx/2,0)+(0,w-dh),EndArrow);
  h-=dh;
  z=(wwx+dd/2-wx/2,h);
  //z=(mx,h);
  filldraw(pic,shift(z)*g,convpen);
  label(pic,"$\{f*g_n\}_{n=0}^{N-1}$",z+(wx/2,w/2),labelpen);

  shipout("cyrc_1exp"+(string) i,pic);
}
