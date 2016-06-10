import graph;
scale(Log);

pen p=linewidth(1);

int mx,my,N;
real k0=1;

struct equil 
{
  real alpha, beta;

  real f(real rho) {
    real sum=0.0;
    for(int i=-mx+1; i < mx; ++i) {
      real i2=i*i;
      for(int j=-my+1; j < my; ++j) {
        real k2=k0*k0*(i2+j*j);
        if(k2 != 0)
          sum += 1.0/(rho+k2);
      }
    }
    return N/sum-rho;
  }

  real fprime(real rho) {
    real sum,sumprime;
    for(int i=-mx+1; i < mx; ++i) {
      real i2=i*i;
      for(int j=-my+1; j < my; ++j) {
        real k2=k0*k0*(i2+j*j);
        if(k2 != 0) {
          real denom=rho+k2;
          sum += 1.0/denom;
          sumprime += 1.0/denom^2;
        }
      }
    }
    return N/sum^2*sumprime-1;
  }

  void operator init(real E, real Z) {
    real r=Z/E;
    write(r);
    
    real rho=1;
    rho=newton(new real(real rho) {return f(rho)-r;},fprime,rho,
               verbose=true);
    beta=0.5*N/(rho*E+Z);
    alpha=rho*beta;
  }
}

int Nx=getint("Nx=");
int Ny=getint("Ny=");

mx=quotient(Nx+1,2);
my=quotient(Ny+1,2);
N=Nx*Ny-1;

file in=input("zkvk").line();
real[][] a=in;
a=transpose(a);

real[] k=a[0];
real[] Zk=a[1];
real[] Ek=Zk/k^2;

draw(graph(k,Ek),p+red);

real E=getreal("E=");
real Z=getreal("Z=");

equil eq=equil(E,Z);

real equiparition(real k) {return pi*k/(eq.alpha+eq.beta*k^2);}
draw(graph(equiparition,min(k),max(k)),p+blue+dashed,
     "$\displaystyle\frac{\pi k}{\alpha+\beta k^2}$");

xaxis("$k$",BottomTop,LeftTicks);
yaxis("$E(k)$",LeftRight,RightTicks);

if(currentpicture.legend.length > 1) attach(legend(),point(E),20E);
