include getparam;
include averages;

import palette;

scale(Linear,Linear);
pen p=linewidth(1);

string Etext="$E$";
string Ztext="$Z$";
string Ptext="$P$";

real[] t,E,Z,P;
real CA=1/2^(1/4);

// Fk=k*fk

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

int k=0;
real G,Lambda;

while(nextrun()) {
  file fin=input(rundir()+"evt").line();
  real[][] a=fin;
  a=transpose(a);
  norm N=fnorm(F);
  real eps=N.f^2;
  real Eta=N.F^2;
  G=sqrt(eps)/nuH^(3/2);                   // Grashof number
  write(G);
  //  Lambda=(N.F/G)^2;             // Lambda := |A^(1/2)f/G|^2
  Lambda=Eta/eps;
  write("Lambda=",Lambda);
  real CG=CA*G;                 // CG := cG in the paper

  real norm=G^2*nuH^2;
  t=a[0]; E=2*a[1]/norm; Z=2*a[2]/norm; P=2*a[3]/norm;

  real Z1=1;
  real Z2=Z1^(1/4)*((3/5)*Z1+(8/5)*Lambda^(1/3)/CG^(4/3))^(3/4);
  real Z3=25*Z2/64;

  int start=getint("start",a[0].length#2,store=false);
  int end=getint("end",a[0].length,store=false);
  t=t[start:end];
  write(t[0],t[t.length-1]);
  E=E[start:end];
  Z=Z[start:end];
  P=P[start:end];
  pen[] p=Rainbow(Z.length);
  
  for(int i=start; i < Z.length; ++i) {
    frame mark;
    fill(mark,scale(0.8mm)*polygon(3+k),p[i-start]);
    add(mark,(Z[i],P[i]));
  }
  
  bool cropx=downcase(getstring("Crop x [Y/n]?","y")) != 'n';
  bool cropy=downcase(getstring("Crop y [Y/n]?","y")) != 'n';

  real Zmin=point(plain.W).x;
  real Zmax=point(plain.E).x;

  if(!cropy) {
  real z1=cropx ? min(Zmax,Z1) : Z1;
  if(Z2 < z1)
    draw(graph(new real(real Z) {return ((4*Lambda*Z1)^(1/3)+1.5*(0.5*CG)^(4/3)*(Z1^(4/3)-Z^(4/3)))^(3/2);},Z2,z1),darkgreen);               //phi1(Z2 to Z1)

  real z2=cropx ? min(Zmax,Z2) : Z2;
  if(Z3 < z2)
  draw(graph(new real(real Z) {return (0.5*CG)^2*(5*(Z*Z2)^(1/2)-4*Z)^2;},Z3,z2),brown);                                          //phi2(Z3 t0 Z2)

  real z3=cropx ? min(Zmax,Z3) : Z3;
  draw(graph(new real(real Z) {return ((2*CG/5)*(6*(Z3^5*Z)^(1/6)-Z))^2;},0,z3),black);                                           //phi3(0 to Z3)
  }

  real Pmax=point(plain.N).y;

  real crop(real x1, real x2=x1) {   
    real bound=1;
    if(cropx) {
      bound=min(bound,x1);
      bound=min(bound,x2);
    }
    return bound;
  }

  //  draw(graph(new real(real Z) {return sqrt(Lambda*Z);},0,Zmax),blue);
  draw(graph(new real(real Z) {return sqrt(Lambda*Z);},0,crop(Pmax^2/Lambda)),blue);
    draw(graph(new real(real Z) {return Z;},0,Zmax),magenta);
  draw(graph(new real(real Z) {return (0.5*CG*Z)^2;},0,crop(2sqrt(Pmax)/CG)),red);
  draw(graph(new real(real Z) {return (2*CG*Z)^2;},0,
             min(Zmax,0.5*sqrt(Pmax)/CG)),pink);
  //             Zmax),pink);



  //draw(graph(new real(real Z) {return kforce^2*Z;},0,min(point(plain.N).y/kforce^2,Zmax)),brown);

  picture bar;
  bounds range=bounds(Zmin,Zmax);
  palette(bar,"$t$",range,(0,0),(0.5cm,6cm),p,NoTicks);
  add(bar.fit(),point(plain.E),30plain.E);

  ++k;
}

xaxis("$2Z$",BottomTop,LeftTicks);
yaxis("$2P$",LeftRight,RightTicks);
