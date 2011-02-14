size(5cm,0);


pair center(path g) {return 0.5*(min(g)+max(g));}



pair b=(1,1), a=(0,0), z=(0,0);
real dh=0.2;
real wx=0.01, wwx=2wx;
real dh=wwx+0.005;
real dd=1.1wwx;
real w=0.75;
real d=wwx+dd, mx=(wwx-wx)/2, h=dh;
pen datapen=yellow, convpen=lightgreen, labelpen=black;
path g, hg, wg;
path p, hp, wp;

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
  size(pic,10cm,0);
  here=0;
  
  datapen=yellow;
  convpen=lightgreen;
  labelpen=black;
  currentpen=black;

  g=box(a,(wx,wx));
  hg=box(a,(wwx,wx));
  wg=box(a,(wwx,wwx));
  pair q1, q2, p1, p2;

  // boxes:
  
  h=dh;
  z=(-1.7wx,h);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"$F$",center(p),labelpen);
  //  p1=(center(p).x,min(p).y);
  p1=center(p)+(wx/2,0);
  
  z=(d+wwx+0.7wx,h);
  p=shift(z)*g;
  filldraw(pic,p,datapen);
  label(pic,"$G$",center(p),labelpen);
  //p2=(center(p).x,min(p).y);
  p2=center(p)-(wx/2,0);
  
  iframe(i);
  //  h-=dh;

  z=(0,h);
  p=shift(z)*g;
  //  q1=(center(p).x,max(p).y);
  q1=center(p)-(wx/2,0);
  draw(pic,p1..q1,EndArrow);
  filldraw(pic,p,datapen);
  label(pic,"$F$",center(p),labelpen);
  wp=shift(z)*wg;
  draw(pic,wp);
  p1=(center(wp).x,min(wp).y);
  
  z=(d,h);
  p=shift(z)*g;
  q2=center(p)+(wx/2,0);
  //q2=(center(p).x,max(p).y);
  draw(pic,p2..q2,EndArrow);
  filldraw(pic,p,datapen);
  label(pic,"$G$",center(p),labelpen);
  wp=shift(z)*wg;
  draw(pic,wp);
  p2=(center(wp).x,min(wp).y);

  /*
  iframe(i);
  h-=dh;
  wg=box(a,(wwx,wwx));
  z=(0,h);
  p=shift(z)*hg;
  q1=(center(p).x,max(p).y);
  draw(pic,p1..q1,EndArrow);
  filldraw(pic,p,datapen);
  wp=shift(z)*wg;
  draw(pic,wp);
  p1=(center(wp).x,min(wp).y);
  
  z=(d,h);
  p=shift(z)*hg;
  q2=(center(p).x,max(p).y);
  draw(pic,p2..q2,EndArrow);
  filldraw(pic,p,datapen);
  wp=shift(z)*wg;
  draw(pic,wp);
  p2=(center(wp).x,min(wp).y);
  */
  
  iframe(i);
  h-=dh;
  z=(0,h);

  wp=shift(z)*wg;
  q1=(center(wp).x,max(wp).y);
  draw(pic,p1..q1,EndArrow);
  filldraw(pic,wp,datapen);
  label(pic,"$f$",center(wp),labelpen);
  p1=(center(wp).x,min(wp).y);
  
  z=(d,h);
  wp=shift(z)*wg;
  q2=(center(wp).x,max(wp).y);
  draw(pic,p2..q2,EndArrow);
  filldraw(pic,wp,datapen);
  label(pic,"$g$",center(wp),labelpen);
  p2=(center(wp).x,min(wp).y);
  
  iframe(i);
  h-=dh;
  //  z=((wwx+dd)/2,h);
  z=(0,h);
  wp=shift(z)*wg;
  q1=(center(wp).x,max(wp).y);
  draw(pic,p1..q1,EndArrow);
  draw(pic,p2..q1,EndArrow);
  filldraw(pic,wp,convpen);
  label(pic,"$fg$",center(wp),labelpen);
  p1=(center(wp).x,min(wp).y);
  
  iframe(i);
  h-=dh;
  //  z=((wwx+dd)/2,h);
  z=(0,h);
  p=shift(z)*g;
  wp=shift(z)*wg;
  q1=(center(wp).x,max(wp).y);
  draw(pic,p1..q1,EndArrow);
  filldraw(pic,p,convpen);
  label(pic,"$F*G$",center(p),labelpen);
  draw(pic,wp);
  //  p1=(center(p).x,min(p).y);
  p1=center(p)-(wx/2,0);

  iframe(i);
  //  h-=dh-(dh-wwx);
  //  z=((wwx+dd)/2-1.7wx,h);
  z-=(1.7wx,0);
  //  z=((wwx+dd)/2,h);
  p=shift(z)*g;
  //  q1=(center(p).x,max(p).y);
  q1=center(p)+(wx/2,0);
  draw(pic,p1..q1,EndArrow);
  filldraw(pic,p,convpen);
  label(pic,"$F*G$",center(p),labelpen);

  
  shipout("conv2psexp"+(string) i,pic);
}
