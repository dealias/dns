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

titlepage("The Fastest Convolution in the West",
	  "John Bowman and Malcolm Roberts","\Blue{University of Alberta}",
	  date="June 16, 2010",
	  url="www.math.ualberta.ca/$\sim$mroberts");

newslide(stepping=false);
title("Outline",newslide=false);
item("Convolution");
subitem("Definition");
subitem("Applications");
item("Fast Convolutions");
subitem("The Convolution Theorem");
//subitem("Cyclic vs.\ Linear");
item("Aliasing Errors");
subitem("Zero-padding");
subitem("Phase-shift dealiasing");
item("Implicit Padding");
subitem("In one dimension");
subitem("In higher dimensions");
item("Centered Hermitian Convolutions");

title("Convolutions");
item("The convolution of the functions $f$ and $g$ is");
equation("(f*g)(t)=\int_{-\infty}^\infty f(\tau) g(t-\tau)\, d\tau.");
item("For example, if $f=g=\chi_{(-1,1)}(t)$");
figure("cyrc_f");
item("Then $f*g$ is:");
figure("cyrc_fg");

title("Applications");
item("Out-of-focus images are a convolution:");
subitem("the actual image is convolved with the aperture opening.");
item("Image filtering:");
subitem("Sobel edge detection is a convolution of the image with a gradient stencil.");
item("Digital signal processing:");
subitem("e.g.\ for low- and high-pass filters.");
item("Correlation analysis.");
item("The Lucas--Lehmer primality test uses fast convolutions.");
subitem("Useful for testing Mersenne primes.");
item("Pseudospectral simulations of fluids:");
subitem("$(u\cdot\del) u$ is a convolution in Fourier space.");

title("Discrete Convolutions");
item("Applications use a {\it discrete linear convolution}:");
equation("{\color{green} (f*g)_n}={\color{green}\sum_{m=0}^n f_m g_{n-m}}.");
item("Calculating $\{(f*g)_n\}_{n=0}^{N-1}$ takes $\O(N^2)$ operations.");
item("The convolution theorem states that convolutions are a multiplications in Fourier space:"); // FIXME: ref?
equation("\mathcal{F}(f*g)=\mathcal{F}(f)\, \mathcal{F}(g)");
remark("where $\{\mathcal{F}(f)\}_k=\sum_{n=0}^{N-1}e^{\frac{2\pi i}{N}kn}f_n$ is the Fourier transform of $\{f_n\}$.");
//remark("where $\mathcal{F}(f)_k=\sum_{n=0}^{N-1}\zeta_N^{nk}f_n$, $\zeta_N=e^{2\pi i/N}$, is the Fourier transform of $\{f_n\}$.");
item("A fast Fourier transform (FFT)  of length $N$ requires $K N \log_2 N$ multiplications~\cite{Gauss1866,Cooley65}.");
item("Convolving using FFTs requires $3KN\log_2 N$ operations.");



title("Cyclic and Linear Convolutions");
item("Fourier transforms map periodic data to periodic data.");
item("Thus, $\mathcal{F}^{-1}[\mathcal{F}(f) \, \mathcal{F}(g) ]$ is a {\it discrete cyclic convolution},");
equation("(f*_{\scriptscriptstyle N}g)_n \doteq \sum_{m=0}^{N-1} f_{m(\mod N)} g_{(n-m) (\mod N)}.");
item("The difference between {\color{green} linear} and cyclic convolutions,");
equation("(f*_{\scriptscriptstyle N}g)_n ={\color{green}\sum_{m=0}^{n} f_m g_{n-m}} + {\color{red} \sum_{m=n+1}^{N-1} f_m g_{n-m+N}},");
step();
remark("is called the {\color{red} \it aliasing error}.");
// FIXME: consider Canuto's dealias.pdf, eq 3.4.9

title("Dealiasing via Explicit Zero-Padding");
item("The cyclic and linear convolutions are equal if we pad $f$ with zeros:");
equation("\{\widetilde f_n\}_{n=0}^{2N-1}=(f_0,f_1,\dots,f_{N-2},f_{N-1},\underbrace{0,\dots,0}_{N})");
// FIXME: include the following?
//item("This can be achieved by choosing $N \ge 2m-1$.");
//item("Since FFT sizes with small prime factors in practice yield the most efficient implementations, the padding is normally extended to $N=2m$.");
item("Then,");
equation("(\widetilde f *_{\scriptscriptstyle 2N}\widetilde g)_n= \sum_{m=0}^{2N-1} \widetilde f_{m (\mod 2N)} \widetilde g_{(n-m) (\mod 2N)},");
step();
equation("= \sum_{m=0}^{N-1} f_{m} \widetilde g_{(n-m) (\mod 2N)},");
step();
equation("= {\color{green} \sum_{m=0}^{n} f_{m} g_{n-m}}.");

title("Dealiasing via Explicit Zero-Padding");
indexedfigure("cyrc_1exp",0,5,"width=22cm");
item("Convolving these padded arrays takes $6K N \log_2 2N$ operations,");
item("and twice the memory of a circular convolution.");
item("CPU speed and memory size have increased much faster than memory bandwidth;  this is the {\it von-Neumann bottleneck}.");
// ref: http://userweb.cs.utexas.edu/~EWD/transcriptions/EWD06xx/EWD692.html

title("Phase-shift Dealiasing");
item("Another possibility is to use a phase shift \cite{Canuto06}.");
item("Define the shifted Fourier transform of $f$ to be");
equation("F^\Delta \doteq \mathcal{F}_k^\Delta(f)=\sum_{n=0}^{N-1} e^{\frac{2\pi i}{N}k(n+\Delta)}f_n.");
item("Then, setting $\Delta=\pi/2$, one has");
equation("f*_{\scriptscriptstyle\Delta}g\doteq {\mathcal{F}^\Delta}^{-1}\left(F^\Delta G^\Delta\right) ={\color{green}\sum_{m=0}^{n} f_m g_{n-m}} - {\color{red} \sum_{m=n+1}^{N-1} f_m g_{n-m+N}},");
remark("which has a dealiasing error with opposite sign:");
step();
//equation("f*g=\frac{1}{2}\left(\mathcal{F}^{-1}\left(F G \right)+{\mathcal{F}^\Delta}^{-1}\left(F^\Delta G^\Delta\right)\right)");
//item("Thus, we can calculate $f*g$ by from two periodic convolutions.");
equation("{\color{green} f*g}=\frac{1}{2}\left(f*_{\scriptscriptstyle N}g +f*_{\scriptscriptstyle\Delta}g\right)");
// note: Canuto's dealias.pdf, eq 3.4.17
item("This requires $6K N \log_2 N$ operations.");
//item("Padding is faster if we need to add fewer than $N$ zeros.");


title("Implicit Padding");
item("Suppose that we want to take a Fourier transform of");
equation("\{f_n\}_{n=0}^{2N-1}, \text{ with }f_n=0 \text{ if } n\geq N");
item("The discrete Fourier transform is a sum:");
equation("\mathcal{F}(f)_k=\sum_{n=0}^{2N-1} e^{\frac{2\pi i}{2N}kn}f_n.");
item("Since $f_n=0$ if $n\geq N$, this is just");
equation("\mathcal{F}(f)_k=\sum_{n=0}^{N-1} e^{\frac{2\pi i}{2N}kn}f_n.");
//equation("\mathcal{F}(f)_k=\sum_{n=0}^{N-1}e^{\frac{ikn}{2N}}f_n.");
item("This is not a FFT, and cannot be done in $\O(N\log_2 N)$.");

title("Implicit Padding");
item("However, if we calculate even and odd terms separately, we get");
equation("\mathcal{F}(f)_{2k}=\sum_{n=0}^{N-1}e^{\frac{2\pi i}{N}kn}f_n, \quad\mathcal{F}(f)_{2k+1}=e^{\frac{ik}{2N}}\sum_{n=0}^{N-1}e^{\frac{2\pi i}{N}kn}f_n,");
step();
remark("which {\it are} FFTs.");
step();
indexedfigure("cyrc_1d",0,2,"width=10cm");
skip();
item("Since Fourier-transformed data is of length $2N$, there are no memory savings.");

title("Implicit Padding");
item("There is one advantage:");
step();
remark("\quad the work buffer is separate from the data buffer.");
step();
indexedfigure("cyrc_1imp",0,3,"width=22cm");
item("The computational complexity is $6 K N log_2 N/2$.");
item("By swapping arrays, we can use out-of-place transforms.");


title("Implicit Padding: speed");
item("The algorithms are comparable in speed:");
figure("ctiming1c","height=13cm");
item("Ours is much more complicated.");

title("Convolutions in Higher Dimensions");
item("An explicitly padded convolution in 2 dimensions requires $12N$ padded FFTs, and 4 times the memory of a cyclic convolution.");
indexedfigure("cyrc_2exp",0,5,"width=14cm");

title("Implicit Convolutions in Higher Dimensions");
item("Implicitly padded 2-dimensional convolutions are done by first doing impliclty padded FFTs in the $x$ direction:");
indexedfigure("cyrc_2dx",0,1,"width=11cm");
item("And then $2N$ one-dimensional convolutions in the $y$-direction:");
indexedfigure("cyrc_2dc",0,1,"width=11cm");


title("Implicit Convolutions in Higher Dimensions");
item("We recover $f*g$ by taking an inverse padded $x$-FFT:");
indexedfigure("cyrc_2dxinv",0,1,"width=6cm");

//newslide();
//item("2D fast convolutions involve a series of FFTs, once for each dimension.");
//item("The first FFT produce needs a separate (but non-contiguous) array:");

//item("$y$-FFTs are done using a 1D work array:");
//figure("cyrc_2dy","height=5cm");
item("An implictly padded convolution in 2 dimensions requires only $9N$ padded FFTs,");
item("and only twice the memory of a cyclic convolution.");
item("The operation count is $6K N \log N/2$.");
//item("Implicit padding uses half the memory of explicit padding in higher dimensions as well.");

//title("Implicit Convolutions in Higher Dimensions");
//item("The transformed arrays are multiplied:");
//figure("cyrc_2dm","height=5cm");
//item("Once we have $F_kG_k$, we take the inverse transform to get $f*g$:");
//figure("cyrc_2dinv","height=4cm");

title("Alternatives");
item("The memory savings could be achieved more simply by using conventional padded transforms.");
item("This requires copying data, which is slow.");
step();
skip();
item("Half of the FFTs in the $x$-direction are on zero-data.");
item("We can skip (``prune\") such transforms:");
figure("cyrc_prune","height=4cm");
item("This is actually slower for large data sets due to memory-striding issues.");

title("Implicit Padding in 2D");
item("Implicit padding is faster in two dimensions:");
figure("ctiming2c","height=13cm");
item("And uses half the memory of explicit padding.");

title("Implicit Padding in 3D");
item("The algorithm is easily extended to three dimensions:");
figure("ctiming3c","height=13cm");
item("Implicit padding uses $1/4$ the memory of explicit padding in 3D.");

title("Centered Hermitian Data");
item("The input $f$ is \emph{centered} if $\{f_n\}_{n=-N/2+1}^{N/2-1} \iff \{F_k\}_{k=-N/2+1}^{N/2-1}$.");
//equation("\{f_n\}_{n=-N/2}^{N/2} \iff \{F_k\}_{k=-N/2}^{N/2}.");
item("If $\{f_n\}$ is real-valued, then $\mathcal{F}(f)$ is \emph{Hermitian}:");
equation("F_{-k}=\conj{F}_k");
//item("Real-to-complex FFTs take $K\frac{N}{2}\log\frac{N}{2}$ multiplies.");
item("The convolution of the centered arrays $f$ and $g$ is");
equation("(f*g)_n = \sum_{p=n-N/2+1}^{N/2-1}f_p g_{n-p}.");
item("Padding centered data increases the array length by $50\%$:");
//figure("cyrc_23","height=4cm"); // FIXME: use this somewhere?
equation("\{\widetilde f_n\}_{n=-N/2+1}^{N-1}
=(f_{-N/2+1},\dots,f_0,\dots,f_{N/2-1},\underbrace{0,\dots,0}_{N/2}).");
item("Phase-shifting is slower than explicit padding for centered data.");

title("Centered Hermitian Data: 1D");
item("The 1D implicit convolution is as fast as explicit padding:");
figure("ctiming1r","height=13cm");
item("And has a comparable memory footprint.");

title("Centered Hermitian Data: 2D");
item("Implicit centered convolutions are faster in higher dimensions:");
figure("ctiming2r","height=13cm");
item("And uses $(2/3)^{d-1}$ the memory in $d$ dimensions.");

//title("Centered Hermitian Data: 3D");
//item("Implicit centered convolutions are faster in higher dimensions:");
//figure("ctiming2r","height=13cm");
// NB: explicit and pruned code not done.

// FIXME:
title("Biconvolutions");
item("The {\it biconvolution\/} of three vectors $F$, $G$, and $H$ is");
equation("\sum_{p=0}^{N-1} \sum_{q=0}^{N-1} F_p G_q H_{k-p-q}.");
item("Computing the transfer function for $Z_4=N^3\sum_\vj \w^4(x_\vj)$ requires
computing the Fourier transform of the cubic quantity~$\w^3$.");
item("This requires a centered Hermitian biconvolution:");
equation("\sum_{p=-m+1}^{m-1}\sum_{q=-m+1}^{m-1}\sum_{r=-m+1}^{m-1}
F_p G_q H_r\d_{p+q+r,k}.");
item("Correctly dealiasing requires a $2/4$ zero padding rule (instead of the usual $2/3$ rule for a single convolution).");

// FIXME
title("2/4 Padding Rule");
item("Computing the transfer function for $Z_4$ with a $2/4$ padding rule means that in a $2048\times 2048$ pseudospectral simulation, the maximum physical wavenumber retained in each direction is only $512$.");
item("For a centered Hermitian biconvolution, implicit padding is twice as fast and uses half of the memory required by conventional explicit padding.");


title("Optimal Problem Sizes");
item("Our main use for this algorithm is pseudo-spectral simulations.");
item("FFTs are faster for highly composite problem sizes:");
subitem("$N=2^n$, $N=3^n$, etc., with $N=2^n$ optimal.");
item("2/3 padding: 341, 683, 1365 etc");
subitem("FFTs have $N=$ 512, 1024, 2048, etc.");
item("Phase-shift dealiasing: $2^n$");
subitem("FFTs are the same size.");
item("Implicit padding: $2^n-1$.");
subitem("sub-transforms are of size $2^{n-1}$.");
item("Implicit padding is optimal for Mersenne-prime sized problem");
//NB: thank Rem for that one


title("Conclusion");
item("Implicitly padded fast convolutions eliminate aliasing errors.");
item("Implicit padding uses $1/2^{d-1}$ the memory of explicit padding in $d$ dimensions.");
item("Computational speedup due to increased data locality and effective pruning.");
//"They use less memory and are faster than explicit zero-padding or phase-shift dealiasing."); // FIXME: check this!
item("Expanding discontiguously easier to program.");
item("``Efficient Dealiased Convolutions without Padding\" submitted to SIAM Journal on Scientific Computing.");
item("A {\tt C++} implementation under the LGPL is available at {\tt http://fftwpp.sourceforge.net/}");
//item("Uses SIMD routines when compiled with the Intel compiler.");
item("Uses the Fastest Fourier Transform in the West ({\tt http://fftw.org/}) for sub-transforms.");
figure("fftw-logo-med.eps","height=1.4cm"); 

//item("Memory savings: in $d$ dimensions implicit padding asymptotically uses $1/2^{d-1}$ of the memory require by conventional explicit padding.");


bibliography("refs");

//  LocalWords:  hypercitebracket dotover cdot mathop nolimits Im vl wl sim
//  LocalWords:  citename Dealiasing dealiasing asyinclude nocite
