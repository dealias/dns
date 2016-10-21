file in;

real[] kb,kc;
real[] F;
int[] Fi,Fj;

void prolog() {
  in=input(run+"/prolog",mode="xdr");
  kb=in.read(1);
  kc=in.read(1);
  while(true) {
    int i=in;
    int j=in;
    real f=in;
    if(eof(in)) break;
    Fi.push(i);
    Fj.push(j);
    F.push(f);
  }
  close(in);
}
