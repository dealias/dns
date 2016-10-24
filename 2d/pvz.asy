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
  G=N.f/nuH^2;                  // Grashof number
  write(G);
  Lambda=(N.F/G)^2;             // Lambda := |A^(1/2)f/G|^2
  real CG=CA*G;                 // CG := cG in the paper

  real norm=nuH^2*G^2;
  write(norm);
  t=a[0]; E=a[1]/norm; Z=a[2]/norm; P=a[3]/norm;

  real Z1=1;
  real Z2=Z1^(1/4)*((3/5)*Z1+(8/5)*Lambda^(1/3)/CG^(4/3))^(3/4);
  real Z3=25*Z2/64;

  int start=getint("start",Z.length#2,store=false);
  start=min(start,Z.length-2);
  pen[] p=BWRainbow(Z.length-start);
  
  for(int i=start; i < E.length; ++i) {
    frame mark;
    fill(mark,scale(0.8mm)*polygon(3+k),p[i-start]);
    add(mark,(E[i],Z[i]));
  }
  
  bool cropx=downcase(getstring("Crop x [Y/n]?","y")) != 'n';
  bool cropy=downcase(getstring("Crop y [Y/n]?","y")) != 'n';

  real crop(real x1, real x2) {   
    real bound=Z1;
    if(cropx) bound=min(bound,x1);
    if(cropy) bound=min(bound,x2);
    return bound;
  }

  //draw(graph(Z,P),yellow,mark);
  real Zmax=point(plain.E).x;
  real Pmax=point(plain.N).y;
  draw(graph(new real(real Z) {return 2*sqrt(Lambda*Z);},0,crop(Zmax,(0.5Pmax)^2/Lambda)),blue); //point(plain.E).x),blue);
  //  draw(graph(new real(real Z) {return Z;},0,Zmax),magenta);
  //draw(graph(new real(real Z) {return (0.5*CG*Z)^2;},0,crop(Zmax,2sqrt(Pmax)/CG)),red);
  //draw(graph(new real(real Z) {return (2*CG*Z)^2;},0,crop(Zmax,0.5*sqrt(Pmax)/CG)),pink);

  if(!cropy) {
  real z1=cropx ? min(Zmax,Z1) : Z1;
  if(Z2 < z1)
    draw(graph(new real(real Z) {return ((4*Lambda*Z1)^(1/3)+1.5*(0.5*CG)^(4/3)*(Z1^(4/3)-Z^(4/3)))^(3/2);},Z2,z1),gray);               //phi1(Z2 to Z1)

  real z2=cropx ? min(Zmax,Z2) : Z2;
  if(Z3 < z2)
  draw(graph(new real(real Z) {return (0.5*CG)^2*(5*(Z*Z2)^(1/2)-4*Z)^2;},Z3,z2),brown);                                          //phi2(Z3 t0 Z2)

  real z3=cropx ? min(Zmax,Z3) : Z3;
  draw(graph(new real(real Z) {return ((2*CG/5)*(6*(Z3^5*Z)^(1/6)-Z))^2;},0,z3),black);                                           //phi3(0 to Z3)
  }
  //draw(graph(new real(real Z) {return kforce^2*Z;},0,min(point(plain.N).y/kforce^2,Zmax)),brown);

  ++k;
}

xaxis("$Z$",BottomTop,LeftTicks);
yaxis("$P$",LeftRight,RightTicks);
