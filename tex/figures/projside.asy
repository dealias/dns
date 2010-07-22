size(12cm,0);

int f=3;
int c=5;
void cross(pair z=(0,0), real w=0.1)
{
  draw((z+w*(1,1))--(z-w*(1,1)));
  draw((z+w*(1,-1))--(z-w*(1,-1)));
}

real d(pair a) {
  return sqrt(a.x^2 + a.y^2);
}

void Arrow(pair a, pair b, real w=0) {
  draw(a-w*(a-b)/d(a-b)..b+w*(a-b)/d(a-b),EndArrow);
}

real h=0.5;
real w=0.05;
real v=1-w;

pair[] fine={(-2,h),(-1,h),(0,h),(1,h)};

dot((0,0),blue);
dot((-2,0),blue);

for (int i=0; i < fine.length; ++i) {
  dot(fine[i],red);
  string lab;
  if (i < 1)
    lab="j-"+(string) (-i+1);
  if (i == 1)
    lab="j";
  if ( i > 1)
    lab="j+"+(string) (i-1);
  label(lab,fine[i]+2*w*N);
}
for (int i=1; i < fine.length; ++i) {
  Arrow(fine[i],(0,0),w);
}

Arrow((-1,h),(-2,0),w);
Arrow((-2,h),(-2,0),w);

Arrow((1,h),(1.5,0.5*h),w);

Arrow((-2.5,0.5*h),(-2,0),w);

label("$1$",(-2,0.3)+w*E);
label("$(1-\beta_j)/2$",0.5*((-1,h)+(-2,0)) +4*w*SE);

label("$1$",0.5*(0,0.6)+w*E);
label("$\beta_{j}/2$",0.5*(-1,h)+4*w*SW);
label("$(1-\beta_{j+2})/2$",0.5*(1,h)+4*w*SE);
real pos=-2.7;
label("fine grid",(pos,h));
label("coarse grid",(pos,0));
