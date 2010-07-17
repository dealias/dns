size(8cm,0);

int n=2;

path g;
real s;
pen p, p2, fillpen, fillpen2, ap;
picture pic;
pair a;
real loffset=4;
real y;
real rh,lh;

int now=0, P;
picture pic;

real dx=2;
picture pic;
string outpre="timestep";
int npics=4;

void go(int thispic) {
  // actually, just have this set all pens to invisible if we're
  // not going to be drawing right now.
  if(thispic > P) {
    p=p2=fillpen=fillpen2=ap=invisible;
    return;
  }
  ap=lightgreen;
  p=black;
  p2=blue;
  fillpen=black;
}


for(P=0; P < npics; ++P){
  go(npics); // set all pens correctly
  write(P);
  pic=new picture;
  size(pic,8cm,0);
  for(int i=0; i < n ; ++i) {
    y=1;
    s=0.1;
    g=scale(s)*unitcircle;
    rh=max(g).x;
    lh=min(g).x;
    a=(dx*i,1);

    if(i>0) go(1);
    filldraw(pic,shift(a)*g,fillpen,p);

    go(0);
    draw(pic,(dx*i+rh,1)..(dx*(i+1)+lh,1),ap,EndArrow);
    
    //second:
    //draw projection
    pair d=(dx*(i-1),0)-a;
    d /= abs(d);
    pair b=(dx*(i-1),0);
    go(2);
    draw(pic,a-d*lh..2*d*lh+b,ap+dashed,EndArrow);
    
    y=0;
    p=blue;
    fillpen2=invisible;
    s *= 2;
    g=scale(s)*unitcircle;
    rh=max(g).x;
    lh=min(g).x;
    
    go(3);
    draw(pic,(dx*i+rh,y)..(dx*(i+1)+lh,y),ap,EndArrow);
    a=(dx*i,0);
    
    //third:
    filldraw(pic,shift(a)*g,fillpen2,p2);
    
    //fourth:
    //draw prolongation
    b=(dx*i,1);
    go(4);
      draw(pic,a+(0,rh)..b+(0,lh/2),ap+dashed,EndArrow);
  }
  shipout(outpre+(string)P,pic);
}
