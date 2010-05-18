orientation=Landscape;

import slide;
import three;

texpreamble("\usepackage[lined,algonl,boxed]{algorithm2e}");

usersetting();

texpreamble("
\def\zhat{\v{\widehat z}}
\let\dotover\dot
\def\dot{{\v \cdot}}
\def\Re{\mathop{\rm Re}\nolimits}
\def\Im{\mathop{\rm Im}\nolimits}
\def\vl{{\v \ell}}
\def\wk{\w_\vk}
\def\wp{\w_\vp}
\def\wq{\w_\vq}
\def\wl{\w_\vl}
\def\Jkpq{\e_{\kpq}}
\def\kpq{{\vk\vp\vq}}
\def\Z{{\mathbb Z}}
\let\ocases\cases

\SetProcFnt{\textnormal}
\setlength{\algomargin}{0.6em}
\SetAlCapSkip{3pt}
\def\ProcNameFnt{\tt}

\def\fft{{\tt fft}}
\def\crfft{{\tt crfft}}
\def\rcfft{{\tt rcfft}}
\def\cconv{{\tt cconv}}
\def\conv{{\tt conv}}
\def\biconv{{\tt biconv}}
\def\build{{\tt build}}
\def\fftpadBackwards{{\tt fftpadBackwards}}
\def\fftpadForwards{{\tt fftpadForwards}}
\def\fftOpadBackwards{{\tt fft0padBackwards}}
\def\fftOpadForwards{{\tt fft0padForwards}}
\def\fftbipadBackwards{{\tt fftbipadBackwards}}
\def\fftbipadForwards{{\tt fftbipadForwards}}
\def\fftObipadBackwards{{\tt fft0bipadBackwards}}
\def\fftObipadForwards{{\tt fft0bipadForwards}}
\SetKwData{xf}{f}
\SetKwData{xu}{u}
\SetKwData{xFk}{F}
\SetKwData{xA}{A}
\SetKwData{xB}{B}
\SetKwData{xg}{g}
\SetKwData{xh}{h}
\SetKwData{xv}{v}
\SetKwData{xw}{w}
\SetKwData{xA}{A}
\SetKwData{xB}{B}
\SetKwData{xC}{C}
\SetKwData{xD}{D}
\SetKwData{xF}{F}
\SetKwData{xG}{G}
\SetKwData{xS}{S}
\SetKwData{xT}{T}
\SetKwData{xU}{U}
\SetKwData{xV}{V}
\SetKwData{xW}{W}
");

usepackage("hypercitebracket");
usepackage("amsmath,mathdef,pde");
texpreamble("\let\cases\ocases");

bibliographystyle("rmp2");

titlepen += darkgreen;

titlepage("The Fastest Convolution in the West",
	  "Malcolm Roberts (University of Alberta)","\Blue{Acknowledgements: John Bowman}",
	  date="May 20, 2010",
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
item("Hermitian Convolutions");

title("Convolutions");
item("The convolution of the functions $f$ and $g$ is");
equation("(f*g)(t)=\int_{-\infty}^\infty f(\tau) g(t-\tau)\, d\tau.");
item("For example, if $f=g=\chi_{(-1,1)}(t)$");
figure("cyrc_f");
item("Then $f*g$ is");
figure("cyrc_fg");

title("Applications");
item("Out-of-focus images are a convolution:");
subitem("the acual image is convolved with the aperture openning.");
item("Image filtering:");
subitem("Sobel edge detection is a convolution of the image with a gradient stencil.");
item("Digital signal processing:");
subitem("e.g.\ for low- and high-pass filters.");
item("Correlation analysis.");
item("The Lucas--Lehmer primality test uses fast convolutions.");
subitem("Useful for testing Mersenne primes.");
item("pseudospectral simulations:");
subitem("$(u\cdot\del) u$ is a convolution in Fourier space.");

title("Discrete Convolutions");
item("Applications typically make use of a discrete convolution:");
equation("(f*g)_n=\sum_{n=p+q} f_p g_q "); //FIXME: choose form
item("Calculationg $(f*g)_n, n=0,\dots N-1$ requires $\O(N^2)$ operations.");
item("The convolution theorem states convolutions are a multitplications in Fourier space:");
equation("\mathcal{F}(f*g)=\mathcal{F}(f)\, \mathcal{F}(g)");
item("Discrete linear convolution sums based on the fast Fourier transform
(FFT) algorithm~\cite{Gauss1866,Cooley65} require $\O(N\log N)$ operations.");


title("Cyclic and Linear Convolutions");
item("FIXME: introduce linear convolutions");
item("Fourier transforms map periodic data to periodic data.");
item("Thus, $\mathcal{F}^{-1}[\mathcal{F}(f) \, \mathcal{F}(g) ]$ is a {\it discrete cyclic convolution}");
equation("\sum_{m=0}^{N-1} f_m g_{n-m},");
remark("where the vectors $f$ and $g$ have period $N$.");
item("The difference between the sum,");
equation("\sum_{m=0}^{n} f_m g_{n-m} -\sum_{m=0}^{N-1} f_m g_{n-m}=
\sum_{m=n+1}^{N-1} f_m g_{n-m}");
remark("are called {\it aliasing errors}.");
// FIXME: consider Canuto's dealias.pdf, eq 3.4.9

title("Dialiasing via Explicit Zero-Pading");
item("If we extend $f_n$ and $g_n$ to be zero for $n\notin (0,\dots N-1)$,
then the cyclic and linear convolution are equal.");
item("This requires $\frac{15}{3} N \log N$ operations,");
item("and $2^d$ times the memory, where $d$ is the dimension.");
item("Memory size and CPU speed have increased much faster than memory bandwidth: the resulting {\it von-Neumann bottleneck} padding a worse choice.");
// FIXME: ref
// http://userweb.cs.utexas.edu/~EWD/transcriptions/EWD06xx/EWD692.html ?


title("Phase-shift Dealiasing");
// FIXME: Canuto's dealias.pdf, eq 3.4.17
item("This requires $15 N \log N$ operations,");
item("but doesn't take up extra memory.");

title("Implict Padding");
item("Suppose that we want to take a Fourier transform of");
equation("f_n, n=0,\dots,2N-1, \text{ with }f_n=0 \text{ if } n\geq N");
item("The discrete Fourier transform is a sum:");
equation("\mathcal{F}(f)_k=\sum_{n=0}^{2N-1}\zeta_{2N}^{kn}f_n.");
item("Since $f_n=0$ if $n\geq N$, this is just");
equation("\mathcal{F}(f)_k=\sum_{n=0}^{N-1}\zeta_{2N}^{kn}f_n.");
item("Unfortunately, this is not an FFT, and cannot be done in $\O(N\log N)$ operations.");

title("Implicit Padding");
item("FIXME");
// so we break it up into even and odd terms
// equation: even and off terms
// These are FFTs:

//equation("u_j\doteq\sum_{k=0}^{N-1}\zeta_N^{jk} U_k\qquad j=0,\ldots,N-1");
// FIXME: the idea is that we can use eq2.1 from dealias.tex
// And it works great
// sub-transforms are done using FFTW
// we can also use out-of-place transforms
// Since the output is twice as big, there are no memory savings.

title("Implict Padding: speed");
item("Unfortunately, there are no speed savings either.");
// FIXME: include timing figure. timing1c.eps
figure("timing1c","height=15cm");

title("Implicit Padding in Higher Dimensions");
item("There is, however, one advantage: the work buffer is separate from the data buffer.");
item("2D fast convolutions involve a series of FFTs, once for each dimension.");
// FIXME: figure
item("The first FFT produces a non-sparse (but non-contiguous) array");
//FIXME: figure
item("Later convolutions may be performed column-by-column.");
// FIXME: figure (movie?) (stepped movie?) (stepped figure!)
item("Doing this with conventional FFTs require copying to a padded buffer.");
// FIXME: figure

title("Implicit Padding in Higher Dimensions");
item("The resulting algorithm has half the memory footprint.");
item("The expected scaling is the same.");
item("However, dut to data locality, it's actually faster:");
// FIXME: 2D timing graph

title("Implicit Padding in Higher Dimensions");
// FIXME: 3D timing graph

title("Hermitian Data");
item("If $f_n=\bar{f}_{N-n}$, the data is {\it Hermitian}"); //FIXME: more like the paper
item("The Fourier and inverser Fourier transforms of Hermitian data are real-valued.");
subitem("As in many applications (e.g.\ pseudo-spectral simulations).");
item("Complex-to-real Fourier transforms scale as $\frac{N}{2}\log\frac{N}{2}$.");

title("Optimal Problem Sizes");
item("Our main use for this algorithm is pseudo-spectral simulations.");
item("FFTs are faster for highly composite problem sizes.");
remark("$N=2^n$, $N=3^n$, etc., with $N=2^n$ optimal.");
item("2/3 padding: 683, etc"); //FIXME: fill this in
remark("FFTs have $N=1024$,etc");
item("Phase-shift dealiasing: $2^n$");
remark("FFTs are the same size");
item("Implicit padding: $2^n-1$.");
remark("sub-transforms are of size $2^{n-1}$");
item("Implicit padding is optimal for Mersenne-prime sized problem");
//NB: thank Rem for that one


title("Conclusion");
//FIXME: yay we're great
//paper in under review
item("Available under the LGPL at:");
remark("{\tt http://fftwpp.sourceforge.net/}");
// uses SIMD routines


if(true){
title("old");


item("Define the $N$th primitive root of unity:");
equation("\zeta_N=\exp\(\fr{2 \pi i}{N}\).");

item("The fast Fourier transform method exploits the properties that
$\zeta_N^r=\zeta_{N/r}$ and $\zeta_N^N=1$.");

item("The unnormalized backwards discrete Fourier transform of\\
$\{f_k: k=0,\ldots,N\}$ is");

equation("f_j\doteq\sum_{k=0}^{N-1}\zeta_N^{jk} F_k\qquad j=0,\ldots,N-1,");

item("The corresponding forward transform is");
equation("F_k\doteq \fr{1}{N}\sum_{j=0}^{N-1}\zeta_N^{-kj} f_j\qquad j=0,\ldots,N-1.");

item("The orthogonality of this transform pair follows from");
equation(
"\sum_{j=0}^{N-1} \zeta_N^{\ell j}=\cases{
N &if $\ell=sN$ for $s\in\Z$,\cr
\ds{\fr{1-\zeta_N^{\ell N}}{1-\zeta_N^{\ell}}}=0 &otherwise.\cr
}");

title("Discrete Linear Convolution");

item("The pseudospectral method requires a {\it linear convolution\/} since wavenumber space is not periodic.");

item("The convolution theorem states:");

equations("\sum_{j=0}^{N-1} f_j g_j\zeta_N^{-jk}
&=&\sum_{j=0}^{N-1} \zeta_N^{-jk}\(\sum_{p=0}^{N-1}\zeta_N^{jp} F_p\)
\(\sum_{q=0}^{N-1}\zeta_N^{jq} G_q\)\endl
&=& \sum_{p=0}^{N-1}\sum_{q=0}^{N-1} F_p G_q \sum_{j=0}^{N-1}\zeta_N^{(-k+p+q)j}\endl
&=& N\sum_s\sum_{p=0}^{N-1} F_p G_{k-p+sN}.");

item("The terms indexed by $s\neq 0$ are called {\it aliases}.");

item("We need to remove the aliases by ensuring that $G_{k-p+sN}=0$
whenever $s\neq 0$.");

item("If $F_p$ and $G_{k-p+sN}$ are nonzero only for $0\le p \le m-1$
and $0\le k-p+sN\le m-1$, then we want $k+sN \le 2m-2$ to have no solutions
for positive $s$.");

item("This can be achieved by choosing $N \ge 2m-1$.");

item("That is, one must {\it zero pad} input data vectors of length~$m$ to length $N\ge 2m-1$.");

item("Physically, {\it explicit zero padding} prevents mode $m-1$ from beating with itself, wrapping around to contaminate mode~$N=0\mod N$");

item("Since FFT sizes with small prime factors in practice yield the most efficient implementations, the padding is normally extended to $N=2m$.");

title("Pruned FFTs");

item("Although explicit padding seems like an obvious waste of memory and computation, the conventional wisdom on avoiding this waste is well summed up
by Steven G.\ Johnson, coauthor of the {\tt FFTW} (``Fastest Fourier Transform in the West\") library \cite{fftwprune}:"); 

remark("\begin{quotation}
{\it
The most common case where people seem to want a pruned FFT is for
zero-padded convolutions, where roughly 50\% of your inputs are zero (to
get a linear convolution from an FFT-based cyclic convolution). Here, a
pruned FFT is \Blue{hardly worth thinking about}, at least in one dimension.
\Green{In higher dimensions, matters change (e.g.\ for a 3d zero-padded array about 1/8 of your inputs are non-zero, and one can fairly easily save a factor of
two or so simply by skipping 1d sub-transforms that are zero).}
}
\end{quotation}");

title("Implicit Padding");

item("If $f_k=0$ for $k \ge m$, one can easily avoid looping over the
unwanted zero Fourier modes by decimating in wavenumber");

equation("
f_{2\ell}
=\ds\sum_{k=0}^{m-1}\zeta_m^{\ell k} F_k,
\qquad
f_{2\ell+1}
=\ds\sum_{k=0}^{m-1}\zeta_m^{\ell k} \zeta_N^kF_k
\qquad
\ell=0,1,\ldots m-1.
");

item("This requires computing two subtransforms, each of size $m$,
for an overall computational scaling of order $2m\log_2 m=N\log_2 m$.");

newslide();

item("Odd and even terms of the convolution can then be computed separately,
multiplied term-by-term, and transformed again to Fourier space:");

equations("
NF_k&=&\sum_{j=0}^{N-1}\zeta_N^{-kj} f_j
=\sum_{\ell=0}^{m-1}\zeta_N^{-k2\ell} f_{2\ell}
+\sum_{\ell=0}^{m-1}\zeta_N^{-k(2\ell+1)} f_{2\ell+1}\endl
&=&\sum_{\ell=0}^{m-1}\zeta_m^{-k\ell} f_{2\ell}
+\zeta_N^{-k}\sum_{\ell=0}^{m-1}\zeta_m^{-k\ell} f_{2\ell+1}
\qquad k=0,\ldots,\fr{N}{2}-1.
");

item("No bit reversal is required at the highest level.");

item("An implicitly padded convolution is implemented as in our {\tt FFTW++} library (version 1.05) as {\tt cconv}({\sf f},{\sf g},{\sf u},{\sf v}) computes an in-place implicitly dealiased convolution of two complex vectors {\sf f} and {\sf g} using two temporary vectors {\sf u} and {\sf v}, each of length~$m$.");

item("This in-place convolution requires six out-of-place transforms, thereby avoiding bit reversal at all levels.");  

if(false) {
remark("
\begin{function}[H]
  \KwIn{vector \xf, vector \xg}
  \KwOut{vector \xf}
  $\xu \leftarrow \fft\inv(\xf)$\;
  $\xv \leftarrow \fft\inv(\xg)$\;
  $\xu \leftarrow \xu * \xv$\;
  \For{$k=0$ \KwTo $m-1$}{
    $\xf[k] \leftarrow \z_{2m}^k\xf[k]$\;
    $\xg[k] \leftarrow \z_{2m}^k\xg[k]$\;
  }
  \medskip
  $\xv \leftarrow \fft\inv(\xf)$\;
  $\xf \leftarrow \fft\inv(\xg)$\;
  $\xv \leftarrow \xv * \xf$\;
  \medskip
  $\xf \leftarrow \fft(\xu)$\;
  $\xu \leftarrow \fft(\xv)$\;
  \medskip
  \For{$k=0$ \KwTo $m-1$}{
    $\xf[k] \leftarrow \xf[k] + \z_{2m}^{-k}\xu[k]$\;
  }
  \Return f/(2m)\;
\end{function}
");
}
title("Implicit Padding in 1D");
figure("timing1c.eps","height=13cm");

title("Implicit Padding in 2D");
figure("timing2c.eps","height=13cm");

title("Implicit Padding in 3D");
figure("timing3c.eps","height=13cm");

title("Hermitian Convolutions");

item("{\it Hermitian convolutions} arise when the input vectors are
Fourier transforms of real data:");

equation("f_{N-k}=\conj{f_k}.");

title("Centered Convolutions");
item("For a {\it centered convolution}, the Fourier origin is at wavenumber zero:");

equation("\sum_{p=k-m+1}^{m-1} f_p g_{k-p}");

item("Here, one needs to pad to $N\ge 3m-2$ to prevent 
mode $m-1$ from beating with itself to contaminate the most negative
(first) mode, corresponding to wavenumber $-m+1$. Since the ratio of the
number of physical to total modes, $(2m-1)/(3m-2)$ is asymptotic to $2/3$
for large $m$, this padding scheme is often referred to as the {\it 2/3
padding rule.}");

item("The Hermiticity condition then appears as");
equation("f_{-k}=\conj{f_k}.");

title("Implicit Hermitician Centered Padding in 1D");
figure("timing1r.eps","height=13cm");

title("Implicit Hermitician Centered Padding in 2D");
figure("timing2r.eps","height=13cm");

title("Biconvolutions");

item("The {\it biconvolution\/} of three vectors $F$, $G$, and $H$ is");

equation("\sum_{p=0}^{N-1} \sum_{q=0}^{N-1} F_p G_q H_{k-p-q}.");

item("Computing the transfer function for $Z_4=N^3\sum_\vj \w^4(x_\vj)$ requires
computing the Fourier transform of the cubic quantity~$\w^3$.");

item("This requires a centered Hermitian biconvolution:");

equation("\sum_{p=-m+1}^{m-1}\sum_{q=-m+1}^{m-1}\sum_{r=-m+1}^{m-1}
F_p G_q H_r\d_{p+q+r,k}.");

item("Correctly dealiasing requires a $2/4$ zero padding rule (instead of the usual $2/3$ rule for a single convolution).");

title("2/4 Padding Rule");
item("Computing the transfer function for $Z_4$ with a $2/4$ padding rule means that in a $2048\times 2048$ pseudospectral simulation, the maximum physical wavenumber retained in each direction is only $512$.");

item("For a centered Hermitian biconvolution, implicit padding is twice as fast and uses half of the memory required by conventional explicit padding.");

title("Implicit Biconvolution in 1D");
figure("timing1b.eps","height=13cm");

title("Implicit Biconvolution in 2D");
figure("timing2b.eps","height=13cm");

title("Conclusions");

item("Memory savings: in $d$ dimensions implicit padding asymptotically uses $1/2^{d-1}$ of the memory require by conventional explicit padding.");

item("Computational savings due to increased data locality: about a factor of two.");

item("Highly optimized versions of these routines have been implemented as a software layer {\tt FFTW++} on top of the {\tt FFTW} library and released under 
the Lesser GNU Public License.");

item("With the advent of this {\tt FFTW++} library, writing a high-performance dealiased pseudospectral code is now a relatively straightforward exercise.");  
}
bibliography("refs");

//  LocalWords:  hypercitebracket dotover cdot mathop nolimits Im vl wl sim
//  LocalWords:  citename Dealiasing dealiasing asyinclude nocite
