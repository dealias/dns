int N=6;

int M=2^N;
int P=13;
assert(4P <= M);

real delta=pi*M/P;

real[] HalfSec=0.5*sequence(new real(int i) {return 1/cos(delta/2^i);},N);
pair[] Exp=sequence(new pair(int i) {return exp(I*delta/2^i);},N);

pair zeta(int k)
{
  int j=N-CTZ(k)-1;
  pair zeta=Exp[j];     
  int I=AND(k,-k);
  Exp[j]=HalfSec[j]*(Exp[j-1]+Exp[N-1-CTZ(2*I+OR(2*I,k-I))]);    
  return zeta;
}
  
write();
write(Exp);
write(HalfSec);

for(int k=1; k < P; ++k) 
  write(zeta(k),expi(2*pi*k/P));


