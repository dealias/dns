size(6cm,0);

int f=3;
int c=5;
real w=0.1;
void cross(pair z=(0,0), real w=0.1)
{
  draw((z+w*(1,1))--(z-w*(1,1)));
  draw((z+w*(1,-1))--(z-w*(1,-1)));
}

pair[] coarse={(-1,0),(0,1),(1,0),(0,-1)};
dot((0,0),red);
label("u",0.1*SE);

for (int i=0; i < coarse.length; ++i) {
  //cross(coarse[i],w);
  dot(coarse[i],blue);
  draw(w*coarse[i]..(1-w)*coarse[i],EndArrow);
  if (i == 0) 
    label("u"+(string) i,coarse[i]+0.1*S);
  else
    label("u"+(string) i,coarse[i]+0.1*E);
  draw(coarse[i]..coarse[(i+1)%4]);
}
label("$\alpha/2$",0.5*coarse[1]+0.1*E);
label("$(1-\alpha)/2$",-0.5*coarse[1]+0.2*E);
label("$\beta/2$",0.5*coarse[0]+0.1*N);
label("$(1-\beta)/2$",-0.5*coarse[0]+0.1*N);


