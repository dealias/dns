file in;

real[] kb,kc;

void prolog() {
  in=input(run+"/prolog",mode="xdr");
  kb=in.read(1);
  kc=in.read(1);
  close(in);
}
