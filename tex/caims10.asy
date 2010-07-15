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
item("Grid geometry:");
subitem("radix-2 decimation");
subitem("radix-4 decimation");
item("Decimated evolution equation");
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
//skip();
item("The large Earth simulator can reach $R=2560$. REF"); // FIXME
item("Under-resolved simulations can have errors at the largest scales.");
indexedfigure("underres_",1,8,"width=16cm");

title("Pseudo-spectral simulations");
item("Spectral simulations evolve the Fourier-transformed Navier--Stokes equations:");
equation("1+1=3.");
item("The advection term $u \del \cdot u$ becomes a convolution
$\sum_{k=p+q} \hat{u}_p (q \cdot \hat{u}_q)$ in Fourier space, taking $\O(n^2)$ operations.");
item("This is best done by transforming back into $x$-space and multiplying, using $\O(n \log n)$ operations.");
item("FIXME: figure of pseudo-spectral"); // FIXME

title("Pseudo-spectral simulations");
item("Aliasing errors can be removed by padding.");
item("Using implicit padding from {\tt fftw++}, this is quite simple:"); // FIXME: ref
item("FIXME: put algorithm for calculating nonlinear source here.");

title("The multispectral method");
item("The large scales are more important but we need the small scales.");
item("We would like to reduce resolution at high wavenumbers.");
item("If $\nu=0$, the system reaches a statistical mechanical equilibrium:
quadratic invariants are evenly distributed between modes.");
equation("E(k)=\alpha, \quad \text{or} \quad E(k)=\frac{1}{\alpha+\beta k^2}.");
item("A variably-decimated grid breaks this equilibrium.");
item("We can use uniformly decimated grids:");
item("FIXME: a figure of decimated grids on top of each other");

title("The multispectral method");
item("The 2D vorticity-based Navier--Stokes equations");
equation("\frac{\partial\omega}{\partial t} + \v{u} \cdot \grad \omega
  =  +\nu \nabla^2 \omega.");
remark("are a good choice on which to develop the method.");
item("The 2D equation is complicated by two invariants, but we need
to evolve just a scalar field, $\omega = \hat{z} \cdot \del \times \v u$.");
item("FIXME: figure of just a plain old grid, showing Hermitian symmetry.");

title("Grid geometry: radix-2");

title("Grid geometry: radix-4");
item("we choose radix-4");

title("Evolution equation on decimated grid");
item("we just change k0");

title("Projecting onto a decimated grid");

title("Prolonging from a decimated grid");

title("time-stepping");
item("explain how we move forward in time");
item("just like in goysr talk.");

title("maybe we have some results? yes? no?");

title("future work");
item("finish coding 2D case");
item("work on 3D case");
item("Develop Runge-Kutta integrators with sub-stage accuracy?");


//bibliography("refs");
