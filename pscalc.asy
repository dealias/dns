// determine nuH and nuL to get steady-state solutions to pseudo-spectral
// simulations.

real C1=3.5;
real C2=3.0;

int pL=getint("Low-wavenumber viscosity degree (power of Laplacian)");
int pH=getint("High-wavenumber viscosity degree (power of Laplacian)",1);


real  kd=0.0, k1=0.0, k2=0.0, eta=1.0;

kd=getreal("kd");
real kf=getreal("forcing wavenumber (kf)");
real deltaf=getreal("range of forcing (deltaf)",1);

eta=getreal("eta");

real gamma=2.0/3.0;
real lambda=gamma;

k1=max(0,kf-0.5*deltaf);
k2=max(0,kf+0.5*deltaf);

real epsilon=eta/kf^2;
real zeta=k2*k2*epsilon/(1.0+(C2/C1)^(1.0/gamma))/kf^2;

real nuH=zeta^(1.0-lambda)*(1.0+(C2/C1)^(1.0/gamma))/(2.0*C2*(
			    (kf^(2.0*pH)-kf^(-4.0/3.0))/(2.0*pH+4.0/3.0)+
			    (kd^(2.0*pH)-kf^(2.0*pH))/(2.0*pH)));


real nuL=(epsilon-zeta/(k2*k2))^(1.0-gamma)*(1.0+(C1/C2)^(1.0/gamma))/
  (2.0*C1*kf^(4.0/3.0)*(
			   (kf^(2.0*pL-2.0)-kf^(-4.0/3.0))
			   /(2.0*pL-2.0/3.0)+
			   (kd^(2.0*pL-2.0)-kf^(2.0*pL-2.0))
			   /(2.0*pL-2.0)));
//real force=sqrt(epsilon/(pi*(k2*k2-k1*k1)));
//write("force="+(string)force);

write("nuH="+(string)nuH);
write("nuL="+(string)nuL);
