import graph;
size(200,IgnoreAspect);

// parameters for the calculation:

// initial conditions with alpha=25, beta=1
int mx=4,my=4;
real E=0.738786091677117;
real Z=5.53034770807208;

int itmax=1000;

real alpha(real beta, real n, real E, real Z)
{
  return (n-beta*Z)/E;
}

real prod(real beta,real[] k2,real E,real Z)
{
  real alpha=alpha(beta,k2.length,E,Z);
  real val=1.0;
  int n=k2.length;
  real alpha=alpha(beta,n,E,Z);
  for(int i=0; i < n; ++i)
    val *= alpha + beta*k2[i];
  return val;
}
real prodskip(real beta, real[] k2, real E, real Z, int s)
{
  real alpha=alpha(beta,k2.length,E,Z);
  real val=1.0;
  int n=k2.length;
  real alpha=alpha(beta,n,E,Z);
  for(int i=0; i < n; ++i)
    if(i != s) val *= alpha + beta*k2[i];
  return val;
}
real prodskip(real beta, real[] k2, real E, real Z, int s, int r)
{
  real alpha=alpha(beta,k2.length,E,Z);
  real val=1.0;
  int n=k2.length;
  real alpha=alpha(beta,n,E,Z);
  for(int i=0; i < n; ++i)
    if(i != s && i != r) val *= alpha + beta*k2[i];
  return val;
}

void setk2(int mx, int my, real[] k2) {
  for(int i=1; i < mx; ++i) {
    int i2=i*i;
    for(int j=0; j < my; ++j) {
      k2.push(i2+j*j);
    }
  }
}

// divide by two since we only look at one quadrant
E *= 0.5;
Z *= 0.5;
real[] k2;
setk2(mx,my,k2);
int n=k2.length;

real f(real beta)
{
  real prod=prod(beta,k2,E,Z);
  real alpha=alpha(beta,n,E,Z);
  real val=E*prod;
  for(int i=0; i < n; ++i) val -= prodskip(beta,k2,E,Z,i);
  return val;
}

real dadb=-Z/E;
real fp(real beta)
{
  real alpha=alpha(beta,n,E,Z);
  real prod=prod(beta,k2,E,Z);
  real ret=0.0;
  for(int i=0; i< n; ++i) 
    ret += (dadb+k2[i])*prodskip(beta,k2,E,Z,i);
  
  ret *= E;
  
  for(int i=0; i < n; ++ i) {
    for(int j=0; j < n; ++ j) {
      if(j!=i) ret -= (dadb+k2[j])*prodskip(beta,k2,E,Z,i,j);
    }
  }
  return ret;
}

{
  picture pic;
  size(pic,200,IgnoreAspect);
  draw(pic,graph(f,0.5,2.1));
  yequals(pic,0,Dotted);
  xaxis(pic,"$\beta$",BottomTop,LeftTicks);
  yaxis(pic,"$f$",LeftRight,RightTicks);
  shipout("root",pic);
}

real beta;
bool bracket=true;
bracket=false;
if(bracket) {
  real lower=-n/E;
  real upper=k2[k2.length-1]*Z/E;
  write(upper);
  beta=newton(itmax,f,fp,lower,upper,true);
} else {
  real betaguess=0.5*n/(E+Z); // correct if alpha=beta
  beta=newton(itmax,f,fp,betaguess,true);
}
write("found beta to be "+(string)beta);

if(beta < 1e100) {
  real alpha=alpha(beta,k2.length,E,Z);
  real e(real alpha,real beta, real ktwo){return 1.0/(alpha+beta*ktwo);}
  real testE=0.0, testZ=0.0;
  for(int i=0; i < k2.length; ++i) {
    testE += e(alpha,beta,k2[i]);
    testZ += k2[i]*e(alpha,beta,k2[i]);
  }
  E *=2.0;
  Z *=2.0;
  write("alpha="+(string)alpha + " beta=" +(string)beta);
  write("original E="+(string) E + ", original Z="+(string) Z);
  write("computed E="+(string) 2testE + ", computed Z="+(string) 2testZ);

  picture pic;
  size(pic,200,IgnoreAspect);
  scale(pic,Log,Log);
  real Ek(real k) {  return 2*pi*k/(alpha+beta*k^2);}
  draw(pic,graph(Ek,sqrt(min(k2)),sqrt(max(k2))));
  xaxis(pic,"$k$",BottomTop,LeftTicks);
  yaxis(pic,"$E(k)$",LeftRight,RightTicks);
  shipout("spec",pic);
}

