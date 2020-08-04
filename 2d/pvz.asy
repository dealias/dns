include getparam;
include averages;

size(400);

import palette;
import ode;

scale(Log,Log);
pen p=linewidth(1);

string Etext="$E$";
string Ztext="$Z$";
string Ptext="$P$";

real[] t,E,Z,P;
real G,Lambda;
int k=0;
real eta,eps,theta;
string tilde;

real CA=1/2^(1/4);

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
  gettime();
  real[][] Tk=transfer();
  eta=2*sum(Tk[ETA]);
  eps=2*sum(Tk[EPS]);
  theta=eta*(N.F/N.f)^2; // Approximate

  string Constant="Constant";
  if(substr(forcing,0,length(Constant)) == Constant) {
    G=N.f/nuH^2;                 // Grashof number for constant forcing
    tilde="";
    Lambda=(N.F/N.f)^2;          // Lambda := |A^(1/2)f|^2/|f|^2
  } else {
    G=sqrt(eps)/nuH^(3/2);       // Grashof number for stochastic forcing
    tilde="\tilde ";
    Lambda=eta/eps;
  }

  write("G=",G);
  write("Lambda=",Lambda);
  real cG=CA*G;

  real norm=G^2*nuH^2;
  t=a[0]; E=2*a[1]/norm; Z=2*a[2]/norm; P=2*a[3]/norm;

  real Z1=1;
  real Z2=Z1^(1/4)*((3/5)*Z1+(4/5)*Lambda^(1/3)/cG^(4/3))^(3/4);
  real Z3=25*Z2/64;

  int start=getint("start",a[0].length#2,store=false);
  int end=getint("end",a[0].length,store=false);
  t=t[start:end];
  real t0=t[0];
  real tmax=t[t.length-1];
  write("t:",t0,tmax);
  E=E[start:end];
  Z=Z[start:end];
  P=P[start:end];
  real incr=(E.length-1)/tmax;
  pen[] p=Rainbow(Z.length);
  
  for(int i=0; i < Z.length; ++i) {
    frame mark;
    fill(mark,scale(0.4mm)*polygon(3+k),p[round(t[i]*incr)]);
    //    fill(mark,scale(0.4mm)*unitcircle,p[round(t[i]*incr)]);
    add(mark,Scale((Z[i],P[i])));
  }
  
  real Zmin=min(Z);
  real Zmax=max(Z);
  real Pmin=min(P);
  real Pmax=max(P);
    real Z1=Z1;
    if(Z2 < Z1)
      draw(graph(new real(real Z) {return ((2*Lambda*Z1)^(1/3)+1.5*(cG/sqrt(2))^(4/3)*(Z1^(4/3)-Z^(4/3)))^(3/2);},Z2,Z1),darkgreen);               //phi1(Z2 to Z1)

    if(Z3 < Z2)
      draw(graph(new real(real Z) {return (cG/sqrt(2))^2*(5*(Z*Z2)^(1/2)-4*Z)^2;},Z3,Z2),brown);                                          //phi2(Z3 t0 Z2)

    draw(graph(new real(real Z) {return (2*sqrt(2)*cG/5*(6*(Z3^5*Z)^(1/6)-Z))^2;},Zmin,Z3),black);                                           //phi3(0 to Z3)

  //  draw(graph(new real(real Z) {return sqrt(Lambda*Z);},0,Zmax),blue);
  //  draw(graph(new real(real Z) {return sqrt(Lambda*Z);},0,Pmax^2/Lambda),blue);
  draw(graph(new real(real Z) {return Z;},Zmin,Zmax),magenta);

  //  draw(graph(new real(real Z) {return (cG*Z)^2/2;},0,sqrt(2*Pmax)/cG),red);
  draw(graph(new real(real Z) {return 8*(cG*Z)^2;},Zmin,sqrt(Pmax/8)/cG),pink);


   draw(graph(new real(real Z) {return (kforce-deltaf/2)^2*Z;},Zmin,min(point(plain.N).y/(kforce-deltaf/2)^2,Zmax)),brown);
   draw(graph(new real(real Z) {return (kforce+deltaf/2)^2*Z;},Zmin,min(point(plain.N).y/(kforce+deltaf/2)^2,Zmax)),brown);

  picture bar;
  bounds range=bounds(Zmin,Zmax);
  palette(bar,"$t$",range,(0,0),(0.5cm,6cm),p,NoTicks);
  add(bar.fit(),point(plain.E),30plain.E);

  ++k;
}

xaxis("$2Z/(\nu "+tilde+"G)^2$",BottomTop,LeftTicks);
yaxis("$2P/(\nu "+tilde+"G)^2$",LeftRight,RightTicks);

