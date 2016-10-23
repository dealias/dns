include getparam;
include averages;

import palette;

scale(Linear,Linear);
pen p=linewidth(1);

string Etext="$E$";
string Ztext="$Z$";
string Ptext="$P$";

real[] t,E,Z,P;
real f,CG,CA=0.5^-0.25;
real G=0,l=0,L=0;
real Z1,Z2,Z3;
int k=0;

real fnorm() {      
  real sum=0;
  for(int k=0; k < F.length; ++k) {
    int i=Fi[k];
    int j=Fj[k];
    real f=F[k];
    l += f;                         // l   := |A^(1/2)f|
    sum += f*f/(i*i+j*j);           // sum := |f|^2    
  }
  return sqrt(sum);                 // |f|
}

void Zi(){
  Z1=1.0;
  Z2=Z1^(1/4)*((3/5.0)*Z1 + ((8/5)*L^(1/6))/(CG^(4/3)))^(3/4);
  Z3=25*Z2/64;
}

while(nextrun()) {
  file fin=input(rundir()+"evt").line();
  real[][] a=fin;
  a=transpose(a);
  G=fnorm()/nuH^2;                  // Grashof number
  write(G);
  L=l/G;                            // L  := Lambda1^(1/2) := |A^(1/2)h|, h=f/G
  CG=CA*G;                          // CG := cG in the paper
  real norm=nuH^2*G^2;              // norm=nuH^2*k0^2/eta;
  write(norm);
  t=a[0]; E=a[1]/norm; Z=a[2]/norm; P=a[3]/norm;
  real E1=a[1][0]/norm;
  //real Lambda=(eta*norm);
  Zi();
  int start=getint("start",Z.length#2,store=false);
  start=min(start,Z.length-2);
  pen[] p=BWRainbow(Z.length-start);
  
  for(int i=start; i < E.length; ++i) {
    frame mark;
    fill(mark,scale(0.8mm)*polygon(3+k),p[i-start]);
    add(mark,(E[i],Z[i]));
  }
  
  //draw(graph(Z,P),yellow,mark);
  draw(graph(new real(real Z) {return 2*L*Z^0.5;},0,1),blue); //point(plain.E).x),blue);
  draw(graph(new real(real Z) {return Z;},0,point(plain.E).x),magenta);
  draw(graph(new real(real Z) {return (0.5*CG*Z)^2;},0,point(plain.E).x),red);
  draw(graph(new real(real Z) {return (2*CG*Z)^2;},0,point(plain.E).x),pink);
  draw(graph(new real(real Z) {return ((4*L^2*Z)^(1/3)+1.5*(0.5*CG)^(4/3)*(Z1^(4/3)-Z^(4/3)))^(3/2);;},Z2,Z1),gray);               //phi1(Z2 to Z1)
  draw(graph(new real(real Z) {return (0.5*CG)^2*(5*(Z*Z2)^(1/2)-4*Z)^2;;},Z3,Z2),brown);                                          //phi2(Z3 t0 Z2)
  draw(graph(new real(real Z) {return ((2*CG/5)*(6*(Z3^5*Z)^(1/6)-Z))^2;;},0,Z3),black);                                           //phi3(0 to Z3)
  //draw(graph(new real(real Z) {return kforce^2*Z;},0,min(point(plain.N).y/kforce^2,point(plain.E).x)),brown);
  ++k;
}

xaxis("$Z$",BottomTop,LeftTicks);
yaxis("$P$",LeftRight,RightTicks);
