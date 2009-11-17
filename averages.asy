import math;
include gettime;
include param;

// Transfer components
int NL=0;
int LIN=1;

real[][] getintegrals(string dir, real T, real Tmax, real exp, int n=1)
{
  real[][] integral;
  int T0=Tindex(T);
  int T1=Tindex(Tmax);
  string prefix=rundir(dir)+"t";
  in=xinput(prefix+(string) T0);
  real[] time=in.read(1);
  real t0=time[0];
  real t1;

  real[][] final;
    
  // final values      
  in=xinput(prefix+(string) T1);
  time=in.read(1);
  t1=time[0];
  for(int i=0; i < n; ++i)
    final[i]=in.read(1);
  close(in);
    
  // values added at steps which are multiples of rezero
  int lastrezero=rezero > 0 ? T1-(T1 % rezero) : 0;
  while(lastrezero > T0) {
    in=xinput(prefix+(string) lastrezero);
    time=in.read(1);
    for(int i=0; i < n; ++i)
      final[i] += (real[]) in.read(1);
    close(in);
    lastrezero -= rezero;
  }
    
  in=xinput(prefix+(string) T0);
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

real[] angleint(real [] A) 
{
  return 2.0*pi*kc*A;
}

real[] Ek;
real[] k;
real[][] M2;

void Ekavg()
{
  write(kb);
  
  M2=getintegrals("ekvk",T,Tmax,2);
  Ek=0.5*angleint(M2[0]);

  k=Ek > 0 ? kc : null;
  Ek=Ek > 0 ? Ek : null;
}

bool nextrun()
{
  ++n;
  int pos=find(runs,",",lastpos);
  if(lastpos == -1) {run=""; return false;}
  run=substr(runs,lastpos,pos-lastpos);
  lastpos=pos > 0 ? pos+1 : -1;
  eval("include \""+rundir()+"param.asy\";",true);

  int N=round(hypot((Nx-1)/2,(Ny-1)/2));
  kc=sequence(1,N);
  kb=kc-0.5;
  kb.push(N+0.5);

  return true;
}