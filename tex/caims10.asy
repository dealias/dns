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

titlepage("The Fastest Convolution in the West",
	  "Malcolm Roberts","\Blue{University of Alberta}",
	  date="2010-07-25",
	  url="www.math.ualberta.ca/$\sim$mroberts");

newslide(stepping=false);
title("Outline",newslide=false);
item("Convolution");
subitem("Definition");
subitem("Applications");
item("Fast Convolutions");
subitem("The Convolution Theorem");
item("Aliasing Errors");
subitem("Zero-padding");
subitem("Phase-shift dealiasing");
subitem("Implicit Padding");
item("Centered Hermitian Convolutions");
item("Ternary Convolutions");

title("Convolutions");
item("The convolution of the functions $f$ and $g$ is");
equation("(f*g)(t)=\int_{-\infty}^\infty f(\tau) g(t-\tau)\, d\tau.");
item("For example, if $f=g=\chi_{(-1,1)}(t)$");
figure("cyrc_f");
item("Then $f*g$ is:");
figure("cyrc_fg");


bibliography("refs");
