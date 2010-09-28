// determine nuH and nuL to get steady-state solutions to pseudo-spectral
// simulations.

real C1=3.5;
real C2=3.0;

int pL=1;
int pH=0;
int pL=getint("Low-wavenumber viscosity degree (power of Laplacian)");
int pH=getint("High-wavenumber viscosity degree (power of Laplacian)");


real  kd=0.0, k1=0.0, k2=0.0, epsilon=1.0;
kd=getreal("kd");
k1=getreal("k1");
k2=getreal("k2");
epsilon=getreal("epsilon");


real gamma=2.0/3.0;
real lambda=gamma;


real zeta=k2*k2*epsilon/(1.0+(C2/C1)^(1.0/gamma));
real kf=0.5*(k1+k2);


real pow(real a, real b) {
  // kludge: original C code used pow, didn't feel like replacing them all
  return a^b;
}

real nuH=zeta^(1.0-lambda)*(1.0+pow(C2/C1,1.0/gamma))/(2.0*C2*(
			    (pow(kf,2.0*pH)-pow(kf,-4.0/3.0))/(2.0*pH+4.0/3.0)+
			    (pow(kd,2.0*pH)-pow(kf,2.0*pH))/(2.0*pH)));

write("nuH="+(string)nuH);

real nuL=pow(epsilon-zeta/(k2*k2),1.0-gamma)*(1.0+pow(C1/C2,1.0/gamma))/
  (2.0*C1*pow(kf,4.0/3.0)*(
			   (pow(kf,2.0*pL-2.0)-pow(kf,-4.0/3.0))
			   /(2.0*pL-2.0/3.0)+
			   (pow(kd,2.0*pL-2.0)-pow(kf,2.0*pL-2.0))
			   /(2.0*pL-2.0)));

write("nuL="+(string)nuL);

real force=sqrt(epsilon/(pi*(k2*k2-k1*k1)));


write("force="+(string)force);
