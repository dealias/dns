// polystrophy: energy, enstrophy, palinstrophy, hyperpalinstrophy

// simpson(new real(real k){return (exp(1/k^2)-1)*k^(-5/3);},1,100000000);
// 0.51285662696

import graph3;
import palette;

pen opacity=opacity(settings.outformat == "html" ? 1.0 : 0.5);

size3(25cm,15cm,15cm,keepAspect=false);

currentprojection=orthographic(dir(72,40));

include getparam;
include averages;

scale(Log,Log,Log);

real[] t,E,Z,P;
real G,Lambda;
real eta,eps,theta;
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
 // Account for reality condition
  N.f=sqrt(2*f2); // |f|
  N.F=sqrt(2*F2); // |A^(1/2)f|
  return N;
}

triple Scale(picture pic=currentpicture, triple v)
{
  return (pic.scale.x.T(v.x),pic.scale.y.T(v.y),pic.scale.z.T(v.z));
}

triple Scale(picture pic=currentpicture, real x, real y, real z) {
  return Scale(pic,(x,y,z));
}

real minE,minZ,minP;
real maxE,maxZ,maxP;

int k=0;
while(nextrun()) {
  file fin=input(rundir()+"evt").line();
  real[][] a=fin;
  a=transpose(a);
  norm N=fnorm(F);
  gettime();
  real[][] Tk=transfer();
  eta=2*sum(Tk[ETA]);
  eps=2*sum(Tk[EPS]);

  G=sqrt(eta*(nuH+nuL))/nuH^2; // Grashof number for stochastic forcing
  tilde="\tilde ";
  Lambda=eta/eps;

  write("G=",G);
  write("Lambda=",Lambda);

  real norm=G^2*nuH^2;
  write("norm=",norm);

  t=a[0]; E=2*a[1]/norm; Z=2*a[2]/norm; P=2*a[3]/norm;

  int start=1;
  int stop=t.length;
  t=t[start:stop];
  E=E[start:stop];
  Z=Z[start:stop];
  P=P[start:stop];

  guide3 g=graph(E,Z,P,operator --);

  int ceilquotient(int a, int b) {return (a+b-1)#b;}

  int N=10;
  pen[] p=Rainbow(ceilquotient(E.length,N));
  triple curr=Scale(E[0],Z[0],P[0]);
  for(int i=0; i < E.length-1; ++i) {
    triple next=Scale(E[i+1],Z[i+1],P[i+1]);
    draw(curr--next,p[i#N]);
    curr=next;
  }

  ++k;

  minE=min(E);
  minZ=min(Z);
  minP=min(P);

  maxE=max(E);
  maxZ=max(Z);
  maxP=max(P);

  // E=Z
  draw(extrude(graph(new triple(real E) {return (E,E,minE);},
                     minE,maxE),
               Scale(minE,minE,maxE^0.25)-
               Scale(minE,minE,minE)),
       blue+opacity);

  // Z=E^0.5
  draw(extrude(graph(new triple(real E) {return (E,sqrt(E),minE);},
                     minE,maxE),
               Scale(minE,sqrt(minE),maxE^0.25)-
               Scale(minE,sqrt(minE),minE)),
       blue+opacity);

 // P=Z
  draw(extrude(graph(new triple(real Z) {return (minE,Z,Z);},
                     minE,sqrt(maxE)),
               Scale(maxE,minE,minE)-
               Scale(minE,minE,minE)),
       red+opacity);

 // P=Z^0.5
  draw(extrude(graph(new triple(real Z) {return (minE,Z,sqrt(Z));},
                     minE,sqrt(maxE)),
               Scale(maxE,minE,sqrt(minE))-
               Scale(minE,minE,sqrt(minE))),
       red+opacity);

 // P=E
  draw(extrude(graph(new triple(real E) {return (E,minE,E);},
                     minE,maxE),
               Scale(minE,sqrt(maxE),minE)-
               Scale(minE,minE,minE)),
       heavygreen+opacity);

 // P=E^0.25
  draw(extrude(graph(new triple(real E) {return (E,minE,E^0.25);},
                     minE,maxE),
               Scale(minE,sqrt(maxE),minE^0.25)-
               Scale(minE,minE,minE^0.25)),
       heavygreen+opacity);

}

xaxis3(Label("$2E/(\nu "+tilde+"G)^2$",position=1),
       YZEquals(minE,minP,extend=false),OutTicks);
yaxis3(Label("$2Z/(\nu "+tilde+"G)^2$",position=1),
       XZEquals(minE,minP,extend=false),minE,sqrt(maxE),OutTicks);
zaxis3(Label("$2P/(\nu "+tilde+"G)^2$",position=1),
       XYEquals(minE,minE,extend=false),OutTicks);
