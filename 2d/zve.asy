include getparam;
include averages;

size(12cm);

import palette;

scale(Log,Log);
defaultpen(linewidth(1));

real[] t,E,Z;
real G;
int k=0;
real eta,eps;
string tilde;

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
  //  write("f2=",f2);
 // Account for reality condition
  N.f=sqrt(2*f2); // |f|
  N.F=sqrt(2*F2); // |A^(1/2)f|
  return N;
}

while(nextrun()) {
  file fin=input(rundir()+"evt").line();
  real[][] a=fin;
  a=transpose(a);
  norm N=fnorm(F);
  gettime();
  real[][] Tk=transfer();
  eta=2*sum(Tk[ETA]);
  eps=2*sum(Tk[EPS]);

  string Constant="Constant";
  if(substr(forcing,0,length(Constant)) == Constant) {
    G=N.f/nuH^2;                 // Grashof number for constant forcing
    tilde="";
  } else {
    G=sqrt(eps*(nuH+nuL))/nuH^2;      // Grashof number for stochastic forcing
    tilde="\tilde ";
  }

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
    add(mark,Scale((E[i],Z[i])));
  }

  real Emin=min(E);
  real Emax=max(E);
  real Zmin=min(Z);
  real Zmax=max(Z);

    real kfm=kforce-deltaf/2;
    real kfp=kforce+deltaf/2;
    write(kfm,kfp);

    draw(graph(new real(real E) {return E;},Emin,1),red);
    //    draw(graph(new real(real E) {return kfp^2*E;},Emin,1/kfp^2),magenta);
    //   draw(graph(new real(real E) {return kfm^2*E;},Emin,1/kfm^2),magenta);
    draw(graph(new real(real E) {return sqrt(E);},Emin,1),blue);

  //  real tau=1-sqrt(2-sqrt(2));
  //  real alpha=tau/sqrt(2);
  //  draw(graph(new real(real E) {return alpha*sqrt(E);},0,crop(Emax,(Zmax/alpha)^2)),0.5*magenta);

  //xequals(nuH^2/4,0.5*green);
  //yequals(nuH^2*kforce^2/4,0.5*blue);

  picture bar;
  bounds range=bounds(Emin,Emax);
  palette(bar,"$t$",range,(0,0),(0.25cm,0.6*currentpicture.ysize),p,NoTicks);
  add(bar.fit(),point(plain.E),15plain.E);

  ++k;
}

xaxis("$2E/(\nu "+tilde+"G)^2$",BottomTop,LeftTicks);
yaxis("$2Z/(\nu "+tilde+"G)^2$",LeftRight,RightTicks);

