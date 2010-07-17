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
int npics=9;

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


for(int M=0; M < 2; ++M) {
  int pf=4;
  if(M==1) {
    outpre="symstep";
    npics=5;
    pf=2;
  }
  
  for(P=0; P < npics; ++P){
    //  go(npics); // set all pens correctly
    write(P);
    pic=new picture;
    size(pic,8cm,0);

    pair b;
    for(int i=0; i < n ; ++i) {
      go(npics); // set all pens correctly
      y=1;
      s=0.1;
      g=scale(s)*unitcircle;
      rh=max(g).x;
      lh=min(g).x;
      a=(dx*i,1);
      
      go(pf*i);
      filldraw(pic,shift(a)*g,fillpen,p);

      go(pf*i+1);
      draw(pic,(dx*i+rh,1)..(dx*(i+1)+lh,1),ap,EndArrow);
      a=(dx*(i+1),1);
      filldraw(pic,shift(a)*g,fillpen,p);

      if(M==1) {
	go(pf*i+2);
	draw(pic,(dx*(i+1),1+lh)..(dx*(i+1),2*rh),ap+dashed,Arrows);
      }
	
      if(M==0) {
	//second:
	//draw projection
	pair d=(dx*i,0)-a;
	d /= abs(d);
	b=(dx*i,0);
	go(pf*i+2);
	draw(pic,a-d*lh..2*d*lh+b,ap+dashed,EndArrow);
      }
	
      y=0;
      s *= 2;
      g=scale(s)*unitcircle;
      rh=max(g).x;
      lh=min(g).x;
      
      if(M==0) {
	//third:
	go(pf*i+3);
	draw(pic,(dx*i+rh,y)..(dx*(i+1)+lh,y),ap,EndArrow);
	a=(dx*(i+1),0);
	filldraw(pic,shift(a)*g,fillpen2,p2);    
      }
      
      a=(dx*i,0);
      go(pf*i);
      filldraw(pic,shift(a)*g,fillpen2,p2);
      if(M==1) {
	go(pf*i+1);
	a=(dx*(i+1),0);
	filldraw(pic,shift(a)*g,fillpen2,p2);
	draw(pic,(dx*i+rh,0)..(dx*(i+1)+lh,0),ap,EndArrow);
      }
	
      if(M==0) {
	//fourth:
	//draw prolongation
	a=(dx*(i+1),0);
	b=(dx*(i+1),1);
	go(pf*i+4);
	draw(pic,a+(0,rh)..b+(0,lh/2),ap+dashed,EndArrow);
      }
    }
    shipout(outpre+(string)P,pic);
  }
}
