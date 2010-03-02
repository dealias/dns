int N=64;

int N4=4N;
//int N4=8N; // Will autorestore table if N is a power of 2.

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
  int i=AND(k,-k);
  if(j == 2) {
    write(HalfSec[j]);
    write(j-1,L-1-CTZ(2*i+OR(2*i,k-i)));
    write(Exp[j-1],Exp[L-1-CTZ(2*i+OR(2*i,k-i))]);
    write(HalfSec[j]*(Exp[j-1]+Exp[L-1-CTZ(2*i+OR(2*i,k-i))]));    
    write();
  }
  if(j == 3) {
    write(HalfSec[j]);
    Exp[j]=i == k ? -I : I;
    write(j-1,L-1-CTZ(2*i+OR(2*i,k-i)));
    write(Exp[j-1],Exp[L-1-CTZ(2*i+OR(2*i,k-i))]);
    write(HalfSec[j]*(Exp[j-1]+Exp[L-1-CTZ(2*i+OR(2*i,k-i))]));    
    write();
    return expi(2*pi*k/N);
  }
  Exp[j]=HalfSec[j]*(Exp[j-1]+Exp[L-1-CTZ(2*i+OR(2*i,k-i))]);    
  return zeta;
}
  
real x=100;
HalfSec[3]=x;
Exp[0]=(1,2.5/x);
Exp[1]=(1,0.5/x);
Exp[2]=(-1,0.5/x);

write(HalfSec);
write();
write(Exp);
write();

for(int k=1; k < N; ++k) 
  write(zeta(k),expi(2*pi*k/N));



