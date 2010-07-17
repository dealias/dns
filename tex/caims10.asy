orientation=Landscape;

import slide;
import three;

texpreamble("\usepackage[lined,algonl,boxed]{algorithm2e}");

usersetting();

usepackage("hypercitebracket");
usepackage("amsmath,mathdef,pde");
texpreamble("\let\cases\ocases");

bibliographystyle("rmp2");

titlepen += darkgreen;

//abstract:
/*
                             The multispectral method
                  Malcolm Roberts, John Bowman, Bruno Eckhardt
                                 University of Alberta
                           Email: mroberts@math.ualberta.ca

Spectral simulations of high Reynolds-number turbulence require a very
large number of active modes, ranging from the largest scale of the
system to the dissipation scale. Since there is no intermediate range
of quiescent modes, numerical techniques cannot rely on a separation
of scales to reduce the stiffness of the system. The multispectral
method generalizes spectral reduction [Bowman, Shadwick, and Morrison,
Phys. Rev. Lett. 83, 5491 (1999)] to evolve time-dependent PDEs on a
hierarchy of decimated grids in Fourier space. The grids can be
decimated so that low-wavenumber modes are preserved. High- Reynolds
number turbulence can then be simulated using far fewer degrees of
freedom than required for a full pseudospectral simulation. In
previous work [Roberts, A Multi- Spectral Decimation Scheme for
Turbulence Simulations, Master’s thesis, University of Alberta
(2006)], the multispectral technique was demonstrated for shell
models. Here, we apply the technique to the two-dimensional
incompressible Navier-Stokes equation, taking care that the projection
and prolongation operators between the grids conserve both energy and
enstrophy. The nonlinear advection term is handled using efficient
algorithms that we have recently developed for computing implicitly
dealiased convolu- tions (http://fftwpp.sourceforge.net).
*/

titlepage("The Multispectral Method",
	  "Bruno Eckhardt, John Bowman, {\bf Malcolm Roberts}",
	  date="2010-07-25",
	  url="www.math.ualberta.ca/$\sim$mroberts");

newslide(stepping=false);
title("Outline",newslide=false);
item("High Reynolds-Number Turbulence");
item("Pseudospectral simulations");
item("Multispectral simulations");
item("Decimated evolution equation");
item("Grid geometry:");
subitem("radix-2 decimation");
subitem("radix-4 decimation");
item("Synchronizing grids");
subitem("Projection");
subitem("Prolongation");
item("Time-stepping");
//item("initial results?");
item("Future work");

title("High Reynolds-Number Turbulence");
//R=2560 from P173Kaneda.pdf (3D), wall-bounded turbulence
// Jupiter: 50m/s = 5e3 cm/s
// Helium: dynamic viscosity: 8.4 μPa·s = 8.4 e-6 g/cm
// density: 0.08988 g/L (HE) or  1.326 g/cm (Jupiter)
// radius: ~ 70 000 km = 7e11 cm
// R ~ 4e+20 ?
item("The incompressible Navier--Stokes equations");
equation(" \ppt{\v{u}} +\v{u}\cdot\grad\v{u} 
  = -\grad P + \nu\nabla^2 \v{u}, \quad \del\cdot \v{u}=0,");
step();
remark("are characterized by the  Reynolds-number $Re$,");
equation("Re=\frac{U L}{\nu}.");
remark("where $U$ and $L$ are characteristic velocity and length scales, and
$\nu$ is the kinematic viscosity.");
item("Energy is dissipated at scales around $\eta_d \sim \nu^{4/3}.$");
item("The number of modes $N$ required ro resolve this grows as");
equation("N \sim Re^{2.25}.");

title("High Reynolds-Number Turbulence");
item("Airplanes: $R \approx 10^6$.");
item("The Earth's atmosphere $R \approx 10^9$.");
//item("Jupiter's red spot has $R \approx 10^{14}$.");
item("Jupiter's atmosphere has $R \approx 10^{20}$.");
item("The large Earth simulator can reach $R=2560$. REF"); // FIXME
item("Under-resolved simulations can have errors at the largest scales.");
indexedfigure("underres_",1,8,"width=16cm");

title("Pseudo-spectral simulations");
item("Spectral simulations evolve the Fourier-transformed Navier--Stokes equations:");
equation("\frac{\partial }{\partial t}\hat{\v{u}}_{\v k}
+ \[(\hat{\v{u}}\cdot \v{k} ) \hat{\v{u}}\]_{\v k}
= \frac{1}{\rho}i \v{k} P - \nu k^2 \hat{\v{u}}_{\v k}");
item("The advection term $u \del \cdot u$ becomes a convolution");
equation("\sum_{\v{k}=\v{p}+\v{q}} \hat{\v{u}}_{\v{p}} (\v{q} \cdot \hat{\v{u}}_{\v{q}})");
remark("in Fourier space, taking $\O(n^2)$ operations.");
item("This is best done by transforming back into $x$-space and multiplying, using $\O(n \log n)$ operations.");
figure("pseudonl","height=2.9cm");
// FIXME: should I show the convolution, as in the equation above?

title("Pseudo-spectral simulations");
item("Aliasing errors can be removed by padding.");
item("Using implicit padding from {\tt fftw++}, this is quite simple:"); // FIXME: ref
item("FIXME: put algorithm for calculating nonlinear source here.");

title("The multispectral method");
item("The 2D vorticity-based Navier--Stokes equations");
equation("\frac{\partial\omega}{\partial t}+ \v{u} \cdot \grad \omega
  =  \nu \nabla^2 \omega.");
remark("are a good choice on which to develop the method.");
item("The 2D equation is complicated by two invariants, but we need
to evolve just a scalar field, $\omega = \hat{z} \cdot \del \times \v u$.");
figure("hermit","width=14cm");

title("The multispectral method");
item("The large scales are more important but we need the small scales.");
item("We would like to decimate at high wavenumbers.");
item("If $\nu=0$, the system reaches a statistical mechanical equilibrium:
quadratic invariants are evenly distributed between modes:");
equation("E(k)=\frac{1}{\alpha+\beta k^2}.");
item("A variably-decimated grid breaks this equilibrium.");
item("We can use uniformly decimated grids:");
figure("figures/grids","width=8cm");

title("Evolution equation on decimated grid"); // FIXME: finish this section
item("we just change k0");
item("small-scale periodicity");
item("reasonable if the high-pass-filtered vorticity field has a
correlcation length less than $1/k0$.");
item("view decimated vortictity as average of undecimated vorticity.");
item("if $k_0 \rightarrow \lambda k_0$, then");
equation("\ppt{\v{u}} +\v{u}\cdot\grad\v{u} 
= -\frac{1}{\rho}\grad P + \nu\nabla^2 \v{u}.
");
item("is sent to");
equation("\ppt{\v{u}} +\v{u}\cdot\lambda\grad\v{u} 
= -\frac{1}{\rho}\lambda\grad P + \nu\lambda^2\nabla^2 \v{u}.
");
item("high-pass filter source term by removeing low-low-low interactions.");

title("Grid geometry: radix-2");
item("Modes with are removed in a checker-board pattern.");
figure("lambdar2rot","width=10cm");
item("The maximum wavenumber is increased by a factor of $\sqrt{2}$.");
item("The overlapping area is a square lattice when symmetry is taken into account.");

title("Grid geometry: radix-4");
item("Modes with are removed along rows and columns.");
figure("lambda2","width=14cm");
item("The maximum wavenumber is increased by a factor of $2$.");
item("The overlap has simple geometry.");
item("High-pass filtering the source term is correspondingly simple.");
item("Synchronizing the grids is also more straightforward.");

title("Synchronizing the grids"); // FIXME:  step?
item("The grids are synchronized via projection and prolongation.");
//item("There are three distinct geometric relationship between modes.");
item("The undecimated and decimated modes might:");
remark("1) be coincident in Fourier-space:");
figure("rad1","height=0.8cm");
remark("2) share the same row/column (only with radix-4 decimation):");
figure("rad4row","height=0.8cm");
remark("3) or the undecimated mode might be between four decimated modes:");
figure("rad4cross","height=3.8cm");

title("Projecting onto a decimated grid"); // FIXME: remove labels from figures
item("Projection set the decimated modes.");
//item("We must populate the modes on the decimated grid."); // FIXME: awkward
//item("FIXME: explain what we do in these three cases.");
item("Case 1 is simple: invariants are automatically conserved.");
figure("rad1","height=0.8cm");
item("For case 2, we need to conserve energy and enstrophy: we have 2
equations and 2 unknowns."); 
figure("rad4row","height=0.8cm");
item("Case 3 is undertermined: 2 equations and 4 unknowns.");
figure("rad4cross","height=3.8cm");
item("One idea is to distribute a linear combination of E and Z.");

title("Prolonging from a decimated grid");
item("Projections sets the undecimated grid.");
item("Case 1 is again simple.");
figure("rad1","height=0.8cm");
item("In case 2, we send some energy up-scale, and some down-scale.");
figure("rad4row","height=0.8cm");
item("We distribute the change from the decimated grid propotionally.");
item("We can deal with case 3 the same way.");
figure("rad4cross","height=3.8cm");

//item("As the undecimated modes have given, so shall they receive.");
// FIXME: flesh this out.
//item("It is important that we modify the vorticity field, and not its source
//term; doing so would ruin equiparition.");

title("Time-stepping");
item("Once we know");
subitem("the decimated evolution equation, and");
subitem("how to synchronize the grids,");
remark("we can move forward in time:");
indexedfigure("timestep",0,4,"height=3.8cm");
//figure("timestep","height=3.8cm");
item("Alternatively, we can synchronize simultaneously:");
figure("timestep","height=3.8cm"); // FIXME
// FIXME: should this be a stepped image?

//title("maybe we have some results? yes? no?");
//item("unfuckinglikely");

title("Conclusions and Future Work");
center("Conclusions:");
item("The Multispectral scheme dramatically reduces the cost of finding solutions to the Navier--Stokes equations.");
item("This technique can be extended to arbitrarily many grids.");
center("Future Work:");
item("finish coding 2D case");
item("accuracy order.");
item("Develop Runge-Kutta integrators with sub-stage accuracy?");
item("develop symmetric projection/prolongation (aka synchronization)");
item("Look into spectral reduction instead of scaling for decimated
evolution equation."); // FIXME: ref
item("work on 3D case");


//bibliography("refs");
