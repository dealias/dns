import graph;
import palette;

size(15cm,IgnoreAspect);

string run=getstring("run");
int N=getint("number of timesteps:",1);
real Dt=getreal("Dt:",1);

scale(Log,Log);
pen p=linewidth(1);

real[] E=new real[N];

import utils;

int n;
real[] k;

real nuH=0.0003;
real nuL=0.2;

real[][] Ek(real[][][] U)
{
  pair[][] u=fft((pair[][]) U[0],1);
  pair[][] v=fft((pair[][]) U[1],1); 

  real norm=1/n^2;
  
  return sequence(new real[](int i) {
      pair[] ui=u[i];
      pair[] vi=v[i];
      real ki=k[i];
      return
        sequence(new real(int j) {
            real k2=ki^2+k[j]^2;
            return (abs(ui[j])^2+abs(vi[j])^2)*norm;},n);},n);
            //          return k2 > 0 ? abs(ki*vi[j]-k[j]*ui[j])^2*norm/k2 : 0;},n);},n);
            //             pair wk=I*(ki*vi[j]-k[j]*ui[j]);
             //             return k2 > 0 ? abs(wk)^2*norm/k2 : 0;},n);},n);
}

for(int t=0; t < N; ++t) {
  real[][][] u=read(Dt*(t+1),"U");

  real e;
  
  n=u[0].length;

  if(t == 0) {
    int h=(n+1)#2;
    k=concat(sequence(0,h-1),sequence(0,h-2)-(h-1));
  }

  real[][] Ek=Ek(u);

  for(int i=0; i < n; ++i) {
    e += sum(Ek[i]);
  }
  E[t]=0.5*e; // TODO: explain (h-1) normalization
  write("["+string(t+1)+"] ",none);
}

write();
write(E);

