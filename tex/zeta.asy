int N=37;

int N4=4N;
int L=CLZ(1)-CLZ(N4);
int M=2^L;
if(M < N4) {++L; M *= 2;}

real delta=pi*M/N;

real[] HalfSec=0.5*sequence(new real(int i) {return 1/cos(delta/2^i);},L);
pair[] Exp=sequence(new pair(int i) {return exp(I*delta/2^i);},L);

pair zeta(int k)
{
  int j=L-CTZ(k)-1;
  pair zeta=Exp[j];     
  int I=AND(k,-k);
  Exp[j]=HalfSec[j]*(Exp[j-1]+Exp[L-1-CTZ(2*I+OR(2*I,k-I))]);    
  return zeta;
}
  
write();
write(Exp);
write(HalfSec);

for(int k=1; k < N; ++k) 
  write(zeta(k),expi(2*pi*k/N));
