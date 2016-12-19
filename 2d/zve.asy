include getparam;
include averages;

import palette;

scale(Linear,Linear);
pen p=linewidth(1);

string Etext="$E$";
string Ztext="$Z$";
string Ptext="$P$";

real[] t,E,Z,P;
real g=0,G=0,l=0,L=0;
int k=0;

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
  G=N.f/nuH^2;                  // Grashof number
  real eps=0.5*N.f^2;
  real Eta=0.5*N.F^2;
  write(eps,Eta/kforce^2,eta/kforce^2);
  G=sqrt(eta/kforce^2)/nuH^(3/2);
  //G=sqrt(eps)/nuH^(3/2);
  write(G);
  //  real norm=0.5*G^2*nuH^2;
  real norm=0.5*G^2*nuH^2;
  t=a[0]; E=a[1]/norm; Z=a[2]/norm;
  int start=getint("start",a[0].length#2,store=false);
  int end=getint("end",a[0].length,store=false);
  t=t[start:end];
  E=E[start:end];
  Z=Z[start:end];
  write(E.length);
  pen[] p=Rainbow(E.length);
  frame mark;
  for(int i=0; i < E.length; ++i) {
    frame mark;
    fill(mark,scale(0.4mm)*polygon(3+k),p[i]);
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

  picture bar;
  bounds range=bounds(Emin,Emax);
  palette(bar,"$t$",range,(0,0),(0.5cm,6cm),p,NoTicks);
  add(bar.fit(),point(plain.E),30plain.E);

  ++k;
}

xaxis("$2E$",BottomTop,LeftTicks);
yaxis("$2Z$",LeftRight,RightTicks);

