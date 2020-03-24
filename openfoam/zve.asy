import graph;
import palette;

size(15cm,IgnoreAspect);

string run=getstring("run");
int N=getint("number of timesteps:",1);

scale(Log,Log);
pen p=linewidth(1);

//bool cropx=downcase(getstring("Crop x [Y/n]?","y")) != 'n';
bool cropx=false;

real crop(real x1, real x2=x1) {
  real bound=1;
  if(cropx) {
    bound=min(bound,x1);
    bound=min(bound,x2);
  }
  return bound;
}

real[] E=new real[N];
real[] Z=new real[N];

real[][][] read(int i, string file) {
  file in=input(run+"/"+string(i)+"/"+file).line();
  while(!eof(in)) {
    string[] s=split(in);
    if(s.length > 0 && s[0] == "internalField") break;
  }
  int n2=in;
  int n=(int) sqrt(n2);
  assert(n*n == n2);

  string s=in;
  assert(s == "(");

  in.line(false);
  triple[][] z=in.dimension(n,n);
  
  real[][] u=sequence(new real[](int i) {
      triple[] zi=z[i];
      return
        sequence(new real(int j) {return zi[j].x;},n);},n);
  real[][] v=sequence(new real[](int i) {
      triple[] zi=z[i];
      return
        sequence(new real(int j) {return zi[j].y;},n);},n);

  u.cyclic=true;
  v.cyclic=true;

  for(int i=0; i < n; ++i)
    u[i].cyclic=v[i].cyclic=true;

  return new real[][][]{u,v};
}

real[][][] f=read(0,"F");

real f2;

real[][] f0=f[0];
real[][] f1=f[1];

pair[][] Fk=fft((pair[][]) f0,1);

int n=f0.length;

/*
for(int i=0; i < n; ++i) {
  for(int j=0; j < n; ++j) {
    if(abs(Fk[i][j]) > 1000)
      write(i,j,abs(Fk[i][j])/n^2);
  }
}
*/

for(int i=0; i < n; ++i) { 
  real[] fx=f0[i];
  real[] fy=f1[i];
  f2 += dot(fx,fx)+dot(fy,fy);
}

real nuH=0.0003;
real nuL=0.2;

real G=sqrt(f2)/nuH^2;
real norm=1/(G^2*nuH^2);

int h=(n+1)#2;
real[] k=concat(sequence(0,h-1),sequence(0,h-2)-(h-1))/n^2;

pair[][] vorticity(real[][][] U)
{
  pair[][] u=fft((pair[][]) U[0],1);
  pair[][] v=fft((pair[][]) U[1],1); 

  u=sequence(new pair[](int i) {
      pair[] ui=u[i];
      pair[] vi=v[i];
      real ki=k[i];
      return
        sequence(new pair(int j) {return I*(ki*vi[j]-k[j]*ui[j]);},n);},n);
  return fft(u,-1);
}

for(int k=0; k < N; ++k) {
  real[][][] u=read(k+1,"U");

  real e,z;
  
  int n=u[0].length;
  pair[][] w=vorticity(u);

  real[][] u0=u[0];
  real[][] u1=u[1];
  for(int i=0; i < n; ++i) {
    real[] ux=u0[i];
    real[] uy=u1[i];
    pair[] wz=w[i];
    e += dot(ux,ux)+dot(uy,uy);
    z += dot(wz,wz).x;
  }
  E[k]=e*norm;
  Z[k]=z*norm;
  write("["+string(k+1)+"] ",none);
}

write();
write(E,Z,sqrt(Z/E));


pen[] P=Rainbow(E.length);

for(int i=0; i < E.length; ++i) {
  frame mark;
  fill(mark,scale(1mm)*unitcircle,P[i]);
  add(mark,Scale((E[i],Z[i])));
}

real kfm=2.23525;
real kfp=3.60475;
  
real Emin=10^point(plain.W).x;
real Emax=10^point(plain.E).x;
real Zmax=10^point(plain.N).y;

draw(graph(new real(real E) {return E;},Emin,crop(Emax)),gray+p);

draw(graph(new real(real E) {return sqrt(E);},Emin,crop(Emax,Zmax^2)),
     blue+p);

draw(graph(new real(real E) {return kfp^2*E;},
           Emin,min(10^point(plain.N).y/kfp^2,10^point(plain.E).x)),magenta);
draw(graph(new real(real E) {return kfm^2*E;},
           Emin,min(10^point(plain.N).y/kfm^2,10^point(plain.E).x)),magenta);


xaxis("$2E/(\nu G)^2$",BottomTop,LeftTicks);
yaxis("$2Z/(\nu G)^2$",LeftRight,RightTicks);
