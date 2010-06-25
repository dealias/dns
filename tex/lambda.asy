size(6cm,0);

int n=4;

//dotfactor=6;

string[] outnames={"lambdar2","lambdar2rot","lambda2"};
real[] lambda={sqrt(2),sqrt(2),2};
pair[] R={(1,0),exp(-pi*I/4),(1,0)};
path g;
picture pic;
pen p;
real L, s;
pair r;
filltype F;


void drawdots() {
  for(int i=-n; i <=n; ++i) {
    for(int j=0; j <=n; ++j) {
      pair a=L*r*(i,j);
      //      dot(a,p,F);
      g=scale(s)*unitcircle;
      filldraw(pic,shift(a)*g,p);
    }
  }
}

for(int i=0; i < lambda.length; ++i) {
  pic = new picture;
  size(pic,10cm,0);
  p=black;
  F=Fill;
  r=1.0;
  L=1.0;
  s=0.1;
  drawdots();

  p=blue+0.01cm;
  F=NoFill;
  L=lambda[i];
  r=R[i];
  s=0.2;
  drawdots();
  shipout(outnames[i],pic);
}

