int N=4;

bool autorestore=false; // Autorestore most of table if N is a power of 2?

int N4;
int offset;
if(autorestore) {
  N4=8*N;
  offset=1;
} else {
  N4=4*N;
  offset=0;
}

int L=CLZ(1)-CLZ(N4);
int M=2^L;
if(M < N4) {++L; M *= 2;}

real delta=pi*M/N;

real[] HalfSec=0.5*sequence(new real(int i) {return 1/cos(delta/2^i);},L);
pair[] Exp=sequence(new pair(int i) {return exp(I*delta/2^i);},L);

write(L);

pair zeta(int k)
{
  int j=L-CTZ(k)-1;
  pair zeta=Exp[j];     
  int i=AND(k,-k);
  Exp[j]=HalfSec[j]*(Exp[j-1]+Exp[L-1-CTZ(2*i+OR(2*i,k-i))]);    
  return zeta;
}
  
real x=-0.5*realEpsilon;

if(N > 3 && M == N4) {
  HalfSec[3+offset]=-0.5/x;
  Exp[0+offset]=(1,5*x);
  Exp[1+offset]=(1,x);
  Exp[2+offset]=(-1,x);
}

write(HalfSec);
write();
write(Exp);
write();

for(int k=1; k < N; ++k) 
  write(zeta(k),expi(2*pi*k/N));
write();

if(autorestore && M == N4) {
  if(N > 3) {
    Exp[0+offset]=(1,5*x);
    Exp[1+offset]=(1,x);
    Exp[2+offset]=(-1,x);
  }
 
  write(Exp);
  write();

  for(int k=1; k < N; ++k) 
    write(zeta(k),expi(2*pi*k/N));
}




