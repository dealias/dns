size(6cm,0);

string[] outnames={"hermit"};
real[] lambda={1};
pair[] R={(1,0)};

int n0=7;
int n;
path g;
picture pic;
pen p, fillpen;
real L, s;
pair r;
filltype F;


void drawdots() {
  for(int i=-n; i <=n; ++i) {
    for(int j= i >= 0 ? 0 : 1 ; j <=n; ++j) {
      pair a=L*r*(i,j);
      g=scale(s)*unitcircle;
      filldraw(pic,shift(a)*g,fillpen,p);
    }
  }
}

for(int i=0; i < lambda.length; ++i) {
  pic = new picture;
  size(pic,10cm,0);
  n=n0;
  p=black;
  fillpen=black;
  F=Fill;
  r=1.0;
  L=1.0;
  s=0.1;
  drawdots();
  real l=n+2;
  draw(pic,(0,0)..(0,l),EndArrow);
  draw(pic,(-l,0)..(l,0),Arrows);
  
  
  shipout(outnames[i],pic);
}

