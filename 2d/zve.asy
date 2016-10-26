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

void fnorm() {      
  real sum=0;
  for(int k=0; k < F.length; ++k) {
    int i=Fi[k];
    int j=Fj[k];
    real f=F[k];
    l += f*f;                         // l   := |A^(1/2)f|
    sum += f*f/(i*i+j*j);           // sum := |f|^2    
  }
  g= sqrt(sum);                     // g   := |f|
}

while(nextrun()) {
  fnorm();
  file fin=input(rundir()+"evt").line();
  real[][] a=fin;
  a=transpose(a);
  G=g/nuH^2;                        // Grashof number
  write(G);
  real norm=nuH^2*G^2;              // norm=nuH^2*k0^2/eta;
  t=a[0]; E=a[1]/norm; Z=a[2]/norm;
  int start=getint("start",E.length#2,store=false);
  start=min(start,E.length-2);
  pen[] p=BWRainbow(E.length-start);
  frame mark;
  for(int i=start; i < E.length; ++i) {
    frame mark;
    fill(mark,scale(0.4mm)*polygon(3+k),p[i-start]);
    add(mark,(E[i],Z[i]));
  }
  //  real kfm=kforce-deltaf/2;
  real kfp=kforce+deltaf/2;
  draw(graph(new real(real E) {return E;},0,1.0), blue); //point(plain.E).x),blue);
  draw(graph(new real(real E) {return kfp^2*E;},
             0,min(point(plain.N).y/kfp^2,point(plain.E).x)),magenta);
  //  draw(graph(new real(real E) {return sqrt(E);},0,(point(plain.N).y)^2),red);

  ++k;
}

xaxis("$E$",BottomTop,LeftTicks);
yaxis("$Z$",LeftRight,RightTicks);
