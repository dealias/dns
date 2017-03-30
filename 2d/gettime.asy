void readt()
{
  t=input(rundir()+"t");
  if(t.length < 2) abort("not enough temporal data");
}

void gettime(bool prompt=true) 
{
  readt();
  string message="T? (<= "+(string) t[t.length-2]+") [";
  if(prompt) {
    Tmax=t[t.length-1];
    T=min(0.5*Tmax,t[t.length-2]);
    pair z;
    do {
    string s=getstring("T",(string) T,message+"%s] ",false);
    z=(s == "-") ? T : (pair) s;
    } while(!initialized(z));
    T=z.x;
    if(z.y != 0) Tmax=z.y;
    T=min(T,t[t.length-2]);
  } else write(message+(string) T+"]");
  int iT=max(find(t <= T,-1),0);
  T=t[iT];
  Tmax=t[max(find(t <= Tmax,-1),1)];
  if(Tmax == T) {--iT; T=t[iT];}
  write("Averaging from T="+(string) T+" to "+(string) Tmax);
}

int Tindex(real T) 
{
  int Tindex=find(t >= T);
  if(Tindex < 0) abort("T="+(string) T+" is too large for data range");
  return Tindex;
}
