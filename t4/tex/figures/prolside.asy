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

pair[] fine={(-1,h),(0,h),(+1,h)};

pair[] coarse={(-1,0),(1,0)};

for (int i=0; i < coarse.length; ++i) {
  dot(coarse[i],blue);
}
label("$U_-$",coarse[0]+2*w*S);
label("$U_+$",coarse[1]+2*w*S);

for (int i=0; i < fine.length; ++i) {
  dot(fine[i],red);
  string lab;
}
label("$u_-$",fine[0]+2*w*N);
label("$u_0$",fine[1]+2*w*N);
label("$u_+$",fine[2]+2*w*N);


for (int i=0; i < coarse.length; ++i) {
  Arrow(coarse[i],fine[2*i],w);
  Arrow(coarse[i],fine[1],w);
}

label("$\alpha_0$",0.5*(-1,0.7)+w*E);
label("$1-\alpha_0$",0.5*(1,0.7)+w*E);


real pos=-1.7;
label("fine grid",(pos,h));
label("coarse grid",(pos,0));
