include getparam;
include averages;

import palette;

scale(Linear,Linear);
pen p=linewidth(1);

real[] t,E,Z;
real g=0,G=0,l=0,L=0;
int k=0;
real eta,eps;
real Eta,Eps;

struct norm {
  real f,F;
}

norm fnorm(real[] F) {      
  real f2,F2;
    
  for(int k=0; k < F.length; ++k) {
    int i=Fi[k];
    int j=Fj[k];
    real Fk=F[k];
    real Fk2=Fk*Fk;
    F2 += Fk2;
    f2 += Fk2/(i*i+j*j);
  }
  norm N;
  N.f=sqrt(f2); // |f|
  N.F=sqrt(F2); // |A^(1/2)f|
  return N;
}


while(nextrun()) {
  file fin=input(rundir()+"evt").line();
  real[][] a=fin;
  a=transpose(a);
  norm N=fnorm(F);
  Eps=N.f^2;
  Eta=N.F^2;
  gettime();
  real[][] Tk=transfer();
  eta=2*sum(Tk[ETA]);
  eps=2*sum(Tk[EPS]);
  G=sqrt(eps)/nuH^(3/2);                   // Grashof number
  write("G=",G);
  real norm=G^2*nuH^2;
  t=a[0]; E=2*a[1]/norm; Z=2*a[2]/norm;
  int start=getint("start",a[0].length#2,store=false);
  int end=getint("end",a[0].length,store=false);
  t=t[start:end];
  real t0=t[0];
  real tmax=t[t.length-1];
  write("t:",t0,tmax);
  E=E[start:end];
  Z=Z[start:end];
  real incr=(E.length-1)/tmax;
  pen[] p=Rainbow(E.length);

  for(int i=0; i < E.length; ++i) {
    frame mark;
    fill(mark,scale(0.4mm)*polygon(3+k),p[round(t[i]*incr)]);
    add(mark,(E[i],Z[i]));
  }

  bool cropx=downcase(getstring("Crop x [Y/n]?","y")) != 'n';
  bool cropy=downcase(getstring("Crop y [Y/n]?","y")) != 'n';

  real Emin=point(plain.W).x;
  real Emax=point(plain.E).x;
  real Zmax=point(plain.N).y;

  real crop(real x1, real x2=x1) {
    real bound=1;
    if(cropx) {
      bound=min(bound,x1);
      bound=min(bound,x2);
    }
    return bound;
  }

    real kfm=kforce-deltaf/2;
    real kfp=kforce+deltaf/2;
  draw(graph(new real(real E) {return E;},0,crop(Emax)),grey); //point(plain.E).x),blue);
     draw(graph(new real(real E) {return kfp^2*E;},
                0,min(point(plain.N).y/kfp^2,point(plain.E).x)),magenta);
     draw(graph(new real(real E) {return kfm^2*E;},
                 0,min(point(plain.N).y/kfm^2,point(plain.E).x)),magenta);
  draw(graph(new real(real E) {return sqrt(E);},0,crop(Emax,Zmax^2)),brown);

  //  real tau=1-sqrt(2-sqrt(2));
  //  real alpha=tau/sqrt(2);
  //  draw(graph(new real(real E) {return alpha*sqrt(E);},0,crop(Emax,(Zmax/alpha)^2)),0.5*magenta);

  //xequals(nuH^2/4,0.5*green);
  //yequals(nuH^2*kforce^2/4,0.5*blue);

  /*
  picture bar;
  bounds range=bounds(Emin,Emax);
  palette(bar,"$t$",range,(0,0),
          (0.6*currentpicture.xsize,0.25cm),Bottom,
          p,NoTicks);

  add(bar.fit(),point(plain.S),30plain.S);
  */

  ++k;
}

xaxis("$2E/(\nu G)^2$",BottomTop,LeftTicks);
yaxis("$2Z/(\nu G)^2$",LeftRight,RightTicks);

