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
subitem("Grid geometry");
subitem("Time-stepping");
subitem("Projection / Prolongation");
item("Test-bed: shell models of turbulence");
item("2D incompressible turbulence");
subitem("Decimation schemes");
subitem("Projection / Prolongation");
item("Future work");

title("High Reynolds-Number Turbulence");
//R=2560 from P173Kaneda.pdf (3D), wall-bounded turbulence
// Jupiter: 50m/s = 5e3 cm/s
// Helium: dynamic viscosity: 8.4 μPa·s = 8.4 e-6 g/cm
// density: 0.08988 g/L (HE) or  1.326 g/cm (Jupiter)
// radius: ~ 70 000 km = 7e11 cm
// R ~ 4e+20 ?
item("K14 theory");
item("We can reach low Reynolds number.");
item("We need to reach high Reynolds numbers.");
item("e.g.\ Jupiter: Re ~ 4e+20");

title("Pseudo-spectral simulations");
item("put pseudo-code for nonlinear source here!");
item("plug fftwpp");

title("The multispectral method");



//bibliography("refs");
