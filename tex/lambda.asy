size(6cm,0);

string[] outnames={"lambda1","lambdar2","lambdar2rot","lambda2"};
real[] lambda={1,sqrt(2),sqrt(2),2};
pair[] R={(1,0),(1,0),exp(-pi*I/4),(1,0)};

int n;
path g;
picture pic;
pen p, fillpen;
real L, s;
pair r;
filltype F;


void drawdots() {
  for(int i=-n; i <=n; ++i) {
    for(int j=0; j <=n; ++j) {
      pair a=L*r*(i,j);
      g=scale(s)*unitcircle;
      filldraw(pic,shift(a)*g,fillpen,p);
    }
  }
}

for(int i=0; i < lambda.length; ++i) {
  pic = new picture;
  size(pic,10cm,0);
  n=4;
  p=black;
  fillpen=black;
  F=Fill;
  r=1.0;
  L=1.0;
  s=0.1;
  drawdots();

  if(i==0) n=8;
  p=blue;
  fillpen=invisible;
  F=NoFill;
  L=lambda[i];
  r=R[i];
  s=0.2;
  drawdots();
  shipout(outnames[i],pic);
}

