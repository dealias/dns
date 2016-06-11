import graph;
size(500);

//scale(Log,Log);

pen p=linewidth(1);

int mx,my,N;
real k0=1;

struct equil 
{
  real alpha, beta;
  real r;

  real f(real rho) {
    real sum;
    for(int i=-mx+1; i < mx; ++i) {
      int i2=i*i;
      for(int j=-my+1; j < my; ++j) {
        int norm2=i2+j*j;
        if(norm2 > 0)
          sum += 1.0/(rho+k0*k0*norm2);
      }
    }
    return N/sum-rho-r;
  }

  real fprime(real rho) {
    real sum,sumprime;
    for(int i=-mx+1; i < mx; ++i) {
      int i2=i*i;
      for(int j=-my+1; j < my; ++j) {
        int norm2=i2+j*j;
        if(norm2 > 0) {
          real factor=1.0/(rho+k0*k0*norm2);
          sum += factor;
          sumprime += factor^2;
        }
      }
    }
    return N/sum^2*sumprime-1;
  }

  void operator init(real E, real Z) {
    r=Z/E;
    
    real rho=newton(f,fprime,1,verbose=true);
    beta=0.5*N/(rho*E+Z);
    alpha=rho*beta;
    write();
    write("alpha=",alpha);
    write("beta=",beta);
  }

  guide graph(real a, real b) {
    return graph(f,a,b);
  }

}

int Nx=getint("Nx=");
int Ny=getint("Ny=");

mx=(Nx+1)#2;
my=(Ny+1)#2;
N=Nx*Ny-1;

real kmax=hypot(mx,my);

file in=input("zkvk").line();
real[][] a=in;
a=transpose(a);

real[] k=a[0];
real[] Zk=a[1];
real[] Ek=Zk/k^2;

real E=getreal("E=");
real Z=getreal("Z=");

equil eq=equil(E,Z);

picture pic;
size(pic,500,IgnoreAspect);
draw(pic,eq.graph(0,kmax));
xaxis(pic,"$\beta$",BottomTop,LeftTicks);
yaxis(pic,"$f$",LeftRight,RightTicks);
shipout(pic);

draw(graph(k,Ek),p+red);

real equiparition(real k) {return pi*k/(eq.alpha+eq.beta*k^2);}
draw(graph(equiparition,min(k),max(k)),p+blue+dashed,
     "$\displaystyle\frac{\pi k}{\alpha+\beta k^2}$");

xaxis("$k$",BottomTop,LeftTicks);
yaxis("$E(k)$",LeftRight,RightTicks);

if(currentpicture.legend.length > 1) attach(legend(),point(E),20E);
