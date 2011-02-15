// TODO:
// add out-of-place convolution diagram?
// add timing graphs!

orientation=Landscape;

import slide;
import three;

texpreamble("
\usepackage[lined,algonl,boxed]{algorithm2e}
\def\F{\mathcal{F}}
\let\ocases\cases
\def\Partial#1#2{\frac{\partial#1}{\partial#2}}
");
//texpreamble("\let\cases\ocases");

usersetting();

usepackage("hypercitebracket");
usepackage("amsmath,mathdef,pde");
texpreamble("\let\cases\ocases");

bibliographystyle("rmp2");

titlepen += darkgreen;

titlepage("The Fastest Convolution in the West",
	  "Malcolm Roberts and John Bowman","\Blue{University of Alberta}",
	  date="2011-02-15",
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
//subitem("In one dimension");
//subitem("In higher dimensions");
item("Centered Hermitian Convolutions");
item("Ternary Convolutions");

title("Convolutions");
item("The convolution of the functions $f$ and $g$ is");
equation("(f*g)(t)=\int_{-\infty}^\infty f(\tau) g(t-\tau)\, d\tau.");
item("For example, if $f=g=\chi_{(-1,1)}(t)$");
figure("figures/cyrc_f");
item("Then $f*g$ is:");
figure("figures/cyrc_fg");

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
skip();
item("Applications use a {\it discrete linear convolution}:");
skip();
equation("{\color{darkgreen} (F*G)_k}={\color{darkgreen}\sum_{\ell=0}^k F_\ell G_{k-\ell}}.");
skip();
item("Calculating $\{(F*G)_k\}_{k=0}^{N-1}$ directly takes $\O(N^2)$ operations.");
skip();
item("Having so many operations introduces significant numerical error.");

title("Fast Convolutions");
//item("The unnormalized backward Fourier transform of $\{F_k\}_{k=0}^{N-1}$ is");
item("The unnormalized backward Fourier transform is");
equation("\{f_n\}_{n=0}^{N-1}=\F^{-1}[F]=\sum_{k=0}^{N-1} \zeta_N^{kn}F_k");
remark("where $\zeta_N=e^{-2\pi/N}$ is the $N^{\rm th}$ root of unity.");
item("The forward transform is");
//item("The forward transform of $\{f_n\}_{n=0}^{N-1}$ is");
equation("\{F_k\}_{k=0}^{N-1}=\F[f]=\frac{1}{N}\sum_{n=0}^{N-1} \zeta_N^{-kn}f_n");
item("This transform relies on the identity");
equation(
"\sum_{j=0}^{N-1} \zeta_N^{\ell j}=
\Large{\{}\begin{matrix} N \text{ if } \ell=sN for s\in\mathbb{Z},
  \\  0 \text{ otherwise.} \end{matrix}");
// FIXME: write cases correctly?

title("Fast Convolutions");
item("Convolutions are multiplications when Fourier-transformed:");
equations("\F[F*G]=\sum_{j=0}^{N-1} f_j g_j\zeta_N^{-jk}
&=&\sum_{j=0}^{N-1} \zeta_N^{-jk}\(\sum_{p=0}^{N-1}\zeta_N^{jp} F_p\)
\(\sum_{q=0}^{N-1}\zeta_N^{jq} G_q\)\endl
&=& \sum_{p=0}^{N-1}\sum_{q=0}^{N-1} F_p G_q \sum_{j=0}^{N-1}\zeta_N^{(-k+p+q)j}\endl
&=& N\sum_s\sum_{p=0}^{N-1} F_p G_{k-p+sN}.");
//item("Convolving using FFTs requires $3KN\log_2 N$ operations.");
item("The terms with $s\neq 0$ are aliases; this is a cyclic convolution:");
equation("\{F *_{\scriptscriptstyle N} G\}_k
= \sum_{\ell =0}^{N-1} F_{\ell \mod N} G_{(k - \ell) \mod N}.");

title("Dealiasing via Explicit Zero-Padding");
item("If we pad $F$ and $G$ with zeroes, we recover the linear convolution:");
equation("\{\widetilde F_k\}_{n=0}^{2N-1}=(F_0,F_1,\dots,F_{N-2},F_{N-1},\underbrace{0,\dots,0}_{N})");
//item("This can be achieved by choosing $N \ge 2m-1$.");
//item("Since FFT sizes with small prime factors in practice yield the most efficient implementations, the padding is normally extended to $N=2m$.");
item("Then,");
equation("(\widetilde F *_{\scriptscriptstyle 2N}\widetilde G)_k=
\sum_{\ell=0}^{2N-1} \widetilde F_{\ell (\mod 2N)} \widetilde G_{(k-\ell) (\mod 2N)},");
step();
equation("= \sum_{\ell=0}^{N-1} F_{\ell} \widetilde G_{(k-\ell) (\mod 2N)},");
step();
equation("= {\color{darkgreen} \sum_{m=0}^{n} F_{\ell} G_{k-\ell}}.");

title("Dealiasing via Explicit Zero-Padding");
indexedfigure("figures/conv1psexp",0,5,"width=22cm");
item("Convolving these padded arrays takes $6K N \log_2 2N$ operations,");
step();
remark("and twice the memory of a circular convolution.");
item("CPU speed and memory size have increased much faster than memory bandwidth;  this is the {\it von-Neumann bottleneck}.");
// ref: http://userweb.cs.utexas.edu/~EWD/transcriptions/EWD06xx/EWD692.html

title("Phase-shift Dealiasing");
item("The shifted Fourier transform \cite{Patterson71} is");
equation("f^\Delta \doteq \{{\F^\Delta}^{-1}[F]\}_k=\sum_{k=0}^{N-1} e^{-\frac{2\pi i}{N}(n+\Delta)k}F_k.");
item("Then, setting $\Delta=1/2$, one has");
equation("(F*_{\scriptscriptstyle\Delta}G)_k
\doteq \F^\Delta\(f^\Delta g^\Delta\)
= {\color{darkgreen}\sum_{\ell=0}^{k} F_\ell G_{k-\ell}}
- {\color{red} \sum_{\ell=k+1}^{N-1} F_\ell G_{k-\ell+N}},");
step();
remark("which has a dealiasing error with opposite sign.");
step();
item("We recover $F*G$ from two periodic convolutions:");
equation("{\color{darkgreen} F*G}=\frac{1}{2}\left(F*_{\scriptscriptstyle N}G +F*_{\scriptscriptstyle\Delta}G\right).");


title("Phase-shift Dealiasing");
indexedfigure("figures/conv1psph",0,4,"width=22cm");
item("We don't need to copy data to a larger buffer first.");
skip();
item("Convolving these padded arrays takes $6K N \log_2 N$ operations,");
skip();
item("The memory footprint is the same as explicit padding.");
skip();
item("Explicit padding is better if we need to add fewer than $N$ zeros.");

title("Implicit Padding");
item("Suppose that we want to take a Fourier transform of");
equation("\{F_k\}_{k=0}^{2N-1}, \text{ with }F_k=0 \text{ if } k\geq N.");
item("The discrete Fourier transform is a sum:");
equation("\F^{-1}(F)_k=\sum_{k=0}^{2N-1} \zeta_{2N}^{kn}F_k.");
item("Since $F_k=0$ if $k\geq N$, this is just");
equation("\F^{-1}(F)_n=\sum_{k=0}^{N-1} \zeta_{2N}^{kn}F_k.");
//equation("\F(f)_k=\sum_{n=0}^{N-1}e^{\frac{ikn}{2N}}f_n.");
item("This is not a Fourier transform: the FFT algorithm doesn't apply.");

title("Implicit Padding");
item("However, if we calculate even and odd terms separately, we get");
equation("
f_{2n}=\sum_{k=0}^{N-1} \zeta_{2N}^{k2n}F_k=\sum_{k=0}^{N-1} \zeta_{N}^{kn}F_k,
");
equation("
f_{2n+1}=\sum_{k=0}^{N-1}\zeta_{2N}^{k(2n+1)}\,F_k
=\sum_{k=0}^{N-1}\zeta_{N}^{kn}\,\(F_k \zeta_{2N}^{k}\),");
step();
remark("which {\it are} Fourier transforms.");
step();
item("The inverse is the sum of two Fourier transforms:");
equation("F_k=\frac{1}{N}\(
\sum_{n=0}^{N-1}\zeta_N^{-kn}f_{2n}+
\zeta_{2N}^k\sum_{k=0}^{N-1}\zeta_N^{-kn}f_{2n+1}
\).");

title("Implicit Padding");
item("Since Fourier-transformed data is of length $2N$, there are no memory savings.");
item("However, the extra memory need not be contiguous: this will be shown to be quite advantageous.");
indexedfigure("figures/conv1psimp",0,3,"width=22cm");
item("The computational complexity is $6 K N \log_2 N/2$.");
//item("By swapping arrays, we can use out-of-place transforms.");
item("The numerical error is similar to explicit padding.");


title("Implicit Padding: speed");
item("The algorithms are comparable in speed:");
figure("timing1c","height=13cm");
item("Ours is much more complicated.");

title("Convolutions in Higher Dimensions");
item("An explicitly padded convolution in 2 dimensions requires $12N$ padded FFTs, and 4 times the memory of a cyclic convolution.");
indexedfigure("figures/conv2psexp",0,5,"width=14cm");

title("Convolutions in Higher Dimensions");
item("Notice that $3/4$ of the transformed arrays are zero.");
skip();
item("It is possible to skip these transforms");
step();
remark("this is called a pruned FFT.");
skip();
item("In the absence of a specially optimized routine for pruned FFTs, it can be faster to simply transform the entire array.");
//figure("figures/cyrc_prune","height=4cm"); // FIXME: include?

title("Implicit Convolutions in Higher Dimensions");
item("One can perform an implicitly-padded 2D convolution by first performing a backward transform in the $x$-direction,");
step();
remark("then performing an implicit 1D convolution in the $y$-direction,");
step();
remark("and then performing a forward transform in the $x$-direction.");
step();
indexedfigure("figures/conv2psimp",0,6,"width=14cm");
item("An implicitly padded convolution in 2 dimensions requires only $9N$ padded FFTs,");
step();
remark("and only twice the memory of a cyclic convolution.");
//item("The operation count is $6K N \log N/2$.");
//item("Implicit padding uses half the memory of explicit padding in higher dimensions as well.");

title("Alternatives");
item("The memory savings could be achieved more simply by using conventional padded transforms.");
skip();
item("This requires copying data, which is slow.");
skip();
item("Phase-shift dealiasing has the same memory footprint as ``1/2\" explicit padding.");


title("Implicit Padding in 2D");
item("Implicit padding is faster in two dimensions:");
figure("timing2c","height=13cm"); // FIXME!
item("And uses half the memory of explicit padding.");

title("Implicit Padding in 3D");
item("The algorithm is easily extended to three dimensions:");
figure("timing3c","height=13cm"); // FIXME!
item("Implicit padding uses $1/4$ the memory of explicit padding in 3D.");

title("Centered Hermitian Data");
item("The input $F$ is \emph{centered} if $\{F_k\}_{k=-N/2+1}^{N/2-1} \iff \{f_n\}_{n=-N/2+1}^{N/2-1}$.");
//equation("\{f_n\}_{n=-N/2}^{N/2} \iff \{F_k\}_{k=-N/2}^{N/2}.");
item("If $\{f_n\}$ is real-valued, then $F$ is \emph{Hermitian}:");
equation("F_{-k}=\conj{F}_k");
//item("Real-to-complex FFTs take $K\frac{N}{2}\log\frac{N}{2}$ multiplies.");
item("The convolution of the centered arrays $f$ and $g$ is");
equation("(F*G)_k = \sum_{\ell=k-N/2+1}^{N/2-1}F_\ell G_{k-\ell}.");
item("Padding centered data use a ``$2/3$\" rule:");
//figure("cyrc_23","height=4cm"); // FIXME: use this somewhere?
equation("\{\widetilde F_k\}_{k=-N/2+1}^{N-1}
=(F_{-N/2+1},\dots,F_0,\dots,F_{N/2-1},\underbrace{0,\dots,0}_{N/2}).");
item("Phase-shifting is slower than explicit padding for centered data.");

title("Centered Hermitian Data: 1D");
item("The 1D implicit convolution is as fast as explicit padding:");
figure("timing1r","height=13cm"); // FIXME!
item("And has a comparable memory footprint.");

title("Centered Hermitian Data: 2D");
item("Implicit centered convolutions are faster in higher dimensions:");
figure("timing2r","height=13cm"); // FIXME!;
item("And uses $(2/3)^{d-1}$ the memory in $d$ dimensions.");

// FIXME:
title("Example: 2D pseudospectral Navier--Stokes");
item("These routines are available in the open-source package {\tt FFTW++}");
item("We need to compute:");
equation("
\Partial{\omega}{t}=-\vu\dot\grad\w
=-(\zhat\cross\grad \del^{-2}\w)\dot \grad \w,
");
remark("which appears in Fourier space as");
equation("
\Partial{\w_\vk}{t} 
=\sum_{\vk=\vp+\vq} \frac{p_xq_y - p_y q_x}{q^2}\w_{\v p}\w_{\v q}.
");
item("The right-hand side of this equation may be computed as");
equation(
"{\tt ImplicitHConvolution2}(i k_x\omega,i k_y\omega,i k_y \omega/k^2,-ik_x\omega/k^2).");

//title("Centered Hermitian Data: 3D");
//item("Implicit centered convolutions are faster in higher dimensions:");
//figure("timing3c","height=13cm"); // FIXME!
// NB: explicit and pruned code not done.

title("Optimal Problem Sizes");
item("We use convolutions in pseudo-spectral simulations:");
//item("The input data is centered and real-valued.");
equation("\partial_t u + u\cdot\del u = -\del P +\nu \del^2 u");
remark("is advanced in Fourier space, with $u\cdot\del u$ calculated in $x$-space.");
item("FFTs are faster for highly composite problem sizes:");
subitem("$N=2^n$, $N=3^n$, etc., with $N=2^n$ optimal.");
item("``$2/3$\" padding: 341, 683, 1365 etc");
subitem("FFTs are of size $N=$ 512, 1024, 2048, etc.");
item("Phase-shift dealiasing: $2^n-1$");
subitem("FFTs are of length $2^{n-1}$.");
//subitem("FFTs are the same size.");
item("Implicit padding: $2^n-1$.");
subitem("sub-transforms are of size $2^{n-1}$.");
//item("Implicit padding is optimal for Mersenne-prime sized problem");

title("Ternary Convolutions");
item("The {\it ternary convolution\/} of three vectors $f$, $g$, and $h$ is");
equation("*\left(F,G,H \right)_{k} =
\sum_{a,b,c\in\{0,\dots,N-1\}}F_a\, G_b\, H_c\,\d_{a+b+c,n}.");
//equation("\left[f,g,h \right]_* =\sum_{m=0}^{N-1} \sum_{p=0}^{N-1} f_m g_p h_{n-m-p}.");
item("Computing the transfer function for $Z_4=N^3\sum_\vj \w^4(x_\vj)$ requires
computing the Fourier transform of $\w^3$.");
item("This requires a centered Hermitian ternary convolution:");
equation("*\left(F,G,H \right)_k =\sum_{a,b,c\in\{-\frac{N}{2}+1,\dots,\frac{N}{2}-1\}}F_a\, G_b\, H_c\,\d_{a+b+c,n}.");
item("Correctly dealiasing requires a ``$2/4$\" padding rule.");
item("Computing $Z_4$ using $2048\times 2048$ pseudospectral modes simulation retains a maximum physical wavenumber of only $512$.");

title("Centered Hermitian Ternary Convolutions: 1D");
item("The 1D implicit ternary convolution is as fast as explicit padding:");
figure("timing1t","height=13cm"); // FIXME!
item("And has a comparable memory footprint.");

title("Centered Hermitian Ternary Convolutions: 2D");
item("Implicit centered ternary convolutions are faster in higher 2D:");
figure("timing2t","height=13cm"); // FIXME!
item("And use $(1/2)^{d-1}$ the memory in $d$ dimensions.");

title("\tt FFTW++");
item("A {\tt C++} implementation, ({\tt FFTW++}, LGPL) is available at {\tt http://fftwpp.sourceforge.net/}.");
item("Fastest Fourier Transform in the West ({\tt http://fftw.org/}) provides sub-transforms.");
item("{\tt FFTw++} will use parallelized sub-transforms when they become available .");
item("Available in {\tt FFTW++}:");
step();
subitem("Non-centered convolutions in 1D, 2D, and 3D,");
step();
subitem("Centered Hermitian convolutions in 1D, 2D, and 3D,");
step();
subitem("Centered Hermitian ternary convolutions in 1D, 2D.");

title("Conclusion");
item("Implicitly padded fast convolutions eliminate aliasing errors.");
skip();
item("Implicit padding uses $(p/q)^{d-1}$ the memory of explicit \mbox{$d$-dimensional} ``$p/q$\" padding.");
skip();
item("Computational speedup from skipping a bit-reversal in the FFT and pruning FFTs efficiently..");
skip();
item("Expanding dis-contiguously is easier to program.");
skip();
item("`Efficient Dealiased Convolutions without Padding\" to appear SIAM Journal on Scientific Computing.");


bibliography("refs");

//  LocalWords:  hypercitebracket dotover cdot mathop nolimits Im vl wl sim
//  LocalWords:  citename Dealiasing dealiasing asyinclude nocite
