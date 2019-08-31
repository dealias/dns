include getparam;
include averages;

size(400);

import palette;
import ode;

scale(Linear,Linear);
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
    add(mark,(Z[i],P[i]));
  }
  
  bool cropx=downcase(getstring("Crop x [Y/n]?","y")) != 'n';
  bool cropy=downcase(getstring("Crop y [Y/n]?","y")) != 'n';

  real Zmin=point(plain.W).x;
  real Zmax=point(plain.E).x;

  if(true) {
    if(!cropy) {
      real z1=cropx ? min(Zmax,Z1) : Z1;
      if(Z2 < z1)
        draw(graph(new real(real Z) {return ((2*Lambda*Z1)^(1/3)+1.5*(cG/sqrt(2))^(4/3)*(Z1^(4/3)-Z^(4/3)))^(3/2);},Z2,z1),darkgreen);               //phi1(Z2 to Z1)

      real z2=cropx ? min(Zmax,Z2) : Z2;
      if(Z3 < z2)
        draw(graph(new real(real Z) {return (cG/sqrt(2))^2*(5*(Z*Z2)^(1/2)-4*Z)^2;},Z3,z2),brown);                                          //phi2(Z3 t0 Z2)

      real z3=cropx ? min(Zmax,Z3) : Z3;
      draw(graph(new real(real Z) {return (2*sqrt(2)*cG/5*(6*(Z3^5*Z)^(1/6)-Z))^2;},0,z3),black);                                           //phi3(0 to Z3)
    }
  }

  real Pmax=point(plain.N).y;

  real crop(real x1, real x2=x1) {   
    real bound=Zmax;
    if(cropx) {
      bound=min(bound,x1);
      bound=min(bound,x2);
    }
    return bound;
  }


  //  draw(graph(new real(real Z) {return sqrt(Lambda*Z);},0,Zmax),blue);
  //  draw(graph(new real(real Z) {return sqrt(Lambda*Z);},0,crop(Pmax^2/Lambda)),blue);
  draw(graph(new real(real Z) {return Z;},0,crop(Zmax)),magenta);

  //  draw(graph(new real(real Z) {return (cG*Z)^2/2;},0,crop(sqrt(2*Pmax)/cG)),red);
  draw(graph(new real(real Z) {return 8*(cG*Z)^2;},0,sqrt(Pmax/8)/cG),pink);


   draw(graph(new real(real Z) {return (kforce-deltaf/2)^2*Z;},0,min(point(plain.N).y/(kforce-deltaf/2)^2,Zmax)),brown);
   draw(graph(new real(real Z) {return (kforce+deltaf/2)^2*Z;},0,min(point(plain.N).y/(kforce+deltaf/2)^2,Zmax)),brown);

  picture bar;
  bounds range=bounds(Zmin,Zmax);
  palette(bar,"$t$",range,(0,0),(0.5cm,6cm),p,NoTicks);
  add(bar.fit(),point(plain.E),30plain.E);


  write("integration test");
  real f1(real Z, real P) {
    Z=-Z;
    //    return 3*(cG/2)^(4/3)*(P*Z)^(1/3);
    return (theta+2*nuL*P+2*(cG*P*Z^(1/4)/sqrt(2))^(4/3))/(2*nuL*Z-P);
  }

  real f2(real Z, real P) {
    Z=-Z;
    //    return 3*(cG/2)^(4/3)*(P*Z)^(1/3);
    real zeta=0;
    return (theta+2*nuL*P-P^2/Z+2*sqrt(2)*(cG*P*Z^(1/4)*zeta)^(4/3))/(2*nuL*Z-P);
  }

  write(integrate(eta,f1,-Z1,-Z2,1e-4,dynamic=true,0.0002,0.0004,RK3BS,verbose=true));
  write(integrate(eta,f2,-Z2,-Z3,1e-4,dynamic=true,0.0002,0.0004,RK3BS,verbose=true));
   write("Pmax=",Pmax);

  ++k;
}

xaxis("$2Z/(\nu "+tilde+"G)^2$",BottomTop,LeftTicks);
yaxis("$2P/(\nu "+tilde+"G)^2$",LeftRight,RightTicks);

