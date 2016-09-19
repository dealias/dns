import math;
include gettime;
include param;
include prolog;

// Spectral components
int EK=0;
int NU=1;

// Transfer components
int NL=0;
int LIN=1;

real[][] getintegrals(string dir, real T, real Tmax, int n=1)
{
  real[][] integral;
  int T0=Tindex(T);
  int T1=Tindex(Tmax);
  string prefix=rundir(dir)+"t";
  in=input(prefix+(string) T0,mode="xdr");
  real[] time=in.read(1);
  real t0=time[0];
  real t1;

  real[][] final;
    
  // final values      
  in=input(prefix+(string) T1,mode="xdr");
  time=in.read(1);
  t1=time[0];
  for(int i=0; i < n; ++i)
    final[i]=in.read(1);
  close(in);
    
  // values added at steps which are multiples of rezero
  int lastrezero=rezero > 0 ? T1-(T1 % rezero) : 0;
  while(lastrezero > T0) {
    in=input(prefix+(string) lastrezero,mode="xdr");
    time=in.read(1);
    for(int i=0; i < n; ++i)
      final[i] += (real[]) in.read(1);
    close(in);
    lastrezero -= rezero;
  }
    
  in=input(prefix+(string) T0,mode="xdr");
  time=in.read(1);
  t0=time[0];
  for(int i=0; i < n; ++i)
    integral[i]=final[i]-(real[]) in.read(1);
  close(in);

  real factor=1.0/(t1-t0);
  for(int i=0; i < n; ++i)
    integral[i] *= factor;
  return integral;
}

real[] Ek;
real[] k,kB;

real[][] moment2()
{
  return getintegrals("ekvk",T,Tmax,2);
}

real[][] transfer() 
{
  return getintegrals("transfer",T,Tmax,2);
}

real[][] transferN() 
{
  return getintegrals("transferN",T,Tmax,2);
}

void Ekavg()
{
  Ek=0.5*moment2()[EK];
  bool[] test=Ek > 0;
  k=test ? kc : null;
  Ek=test ? Ek : null;
  test.insert(0,true);
  kB=test ? kb : null;
}

bool nextrun()
{
  ++n;
  int pos=find(runs,",",lastpos);
  if(lastpos == -1) {run=""; return false;}
  run=substr(runs,lastpos,pos-lastpos);
  lastpos=pos > 0 ? pos+1 : -1;
  eval("include \""+rundir()+"param.asy\";",true);
  
  prolog();
  k=kc;
  kB=kb;

  return true;
}

real kd(real fraction=0.5) {
  real[] eta=2*0.5*moment2()[NU];
  real[] ratio=partialsum(eta)/sum(eta);
  return exp(interpolate(ratio,log(kb),fraction));
}
