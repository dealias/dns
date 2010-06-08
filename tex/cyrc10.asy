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
item("Centered Hermitian Convolutions");

title("Convolutions");
item("The convolution of the functions $f$ and $g$ is");
equation("(f*g)(t)=\int_{-\infty}^\infty f(\tau) g(t-\tau)\, d\tau.");
item("For example, if $f=g=\chi_{(-1,1)}(t)$");
figure("cyrc_f");
item("Then $f*g$ is");
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
equation("{\color{green} (f*g)_n}={\color{green}\sum_{m=0}^n f_m g_{n-m}}");
item("Calculating $\{(f*g)_n\}_{n=0}^{N-1}$ takes $\O(N^2)$ operations.");
item("The convolution theorem states that convolutions are a multiplications in Fourier space:");
equation("\mathcal{F}(f*g)=\mathcal{F}(f)\, \mathcal{F}(g)");
remark("where $\{\mathcal{F}(f)\}_k=\sum_{n=0}^{N-1}e^{\frac{2\pi i}{N}kn}f_n$ is the Fourier transform of $\{f_n\}$.");
//remark("where $\mathcal{F}(f)_k=\sum_{n=0}^{N-1}\zeta_N^{nk}f_n$, $\zeta_N=e^{2\pi i/N}$, is the Fourier transform of $\{f_n\}$.");
item("A fast Fourier transform (FFT)  of length $N$ requires $K N \log_2 N$ multiplications~\cite{Gauss1866,Cooley65}.");
item("Convolving using FFTs requires $3KN\log_2 N$ operations.");



title("Cyclic and Linear Convolutions");
item("Fourier transforms map periodic data to periodic data.");
item("Thus, $\mathcal{F}^{-1}[\mathcal{F}(f) \, \mathcal{F}(g) ]$ is a {\it discrete cyclic convolution},");
equation("(f*_{\scriptscriptstyle N}g)_n \doteq \sum_{m=0}^{N-1} f_{m(\mod N)} g_{(n-m) (\mod N)}.");
//remark("where the vectors $f$ and $g$ have period $N$.");
item("The difference between {\color{green} linear} and cyclic convolutions,");
equation("(f*_{\scriptscriptstyle N}g)_n ={\color{green}\sum_{m=0}^{n} f_m g_{n-m}} + {\color{red} \sum_{m=n+1}^{N-1} f_m g_{n-m+N}},");
step();
remark("is called the {\color{red} \it aliasing error}.");
// FIXME: consider Canuto's dealias.pdf, eq 3.4.9

title("Dealiasing via Explicit Zero-Padding");
item("The cyclic and linear convolutions are equal if we pad $f$ with zeros:");
equation("\{\widetilde f_n\}_{n=0}^{2N-1}=(f_0,f_1,\dots,f_{N-2},f_{N-1},\underbrace{0,\dots,0}_{N})");
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
//item("Explicit zero-padding seems wasteful.");
// ref: http://userweb.cs.utexas.edu/~EWD/transcriptions/EWD06xx/EWD692.html

title("Phase-shift Dealiasing");
item("Another possibility is to use a phase shift \cite{Canuto06}.");
item("Define the shifted Fourier transform of $f$ to be");
equation("F^\Delta \doteq \mathcal{F}_k^\Delta(f)=\sum_{n=0}^{N-1} e^{\frac{2\pi i}{N}k(n+\Delta)}f_n,");
item("Then, setting $\Delta=\pi/2$, one has");
equation("f*_{\scriptscriptstyle\Delta}g\doteq {\mathcal{F}^\Delta}^{-1}\left(F^\Delta G^\Delta\right) ={\color{green}\sum_{m=0}^{n} f_m g_{n-m}} - {\color{red} \sum_{m=n+1}^{N-1} f_m g_{n-m+N}}.");
remark("which has a dealiasing error with opposite sign.");
//equation("f*g=\frac{1}{2}\left(\mathcal{F}^{-1}\left(F G \right)+{\mathcal{F}^\Delta}^{-1}\left(F^\Delta G^\Delta\right)\right)");
item("Thus, we can calculate $f*g$ by from two periodic convolutions.");
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
figure("timing1c","height=13cm"); // TODO: redo (with N instead of m?)
item("Ours is much more complicated.");

title("Convolutions in Higher Dimensions");
item("An explicitly padded convolution in 2 dimensions requires $12N$ padded FFTs, and 4 times the memory of a cyclic convolution.");
indexedfigure("cyrc_2exp",0,5,"width=14cm");

title("Implicit Convolutions in Higher Dimensions");
item("Implicitly padded 2-dimensional convolutions are done by first doing impliclty padded FFTs in the $x$ direction:");
indexedfigure("cyrc_2dx",0,1,"width=11cm");
item("And then $2N$ one-dimensional convolutions in the $y$-direction:");
step();
indexedfigure("cyrc_2dc",0,1,"width=11cm");


title("Implicit Convolutions in Higher Dimensions");
item("We recover $f*g$ by taking an inverse padded $x$-FFT:");
//figure("cyrc_2dxinv","width=6cm"); // FIXME: indexedfigure
indexedfigure("cyrc_2dxinv",0,1,"width=6cm");

//newslide();
//item("2D fast convolutions involve a series of FFTs, once for each dimension.");
//item("The first FFT produce needs a separate (but non-contiguous) array:");

//item("$y$-FFTs are done using a 1D work array:");
//figure("cyrc_2dy","height=5cm");
item("An implictly padded convolution in 2 dimensions requires $9N$ padded FFTs, and twice the memory of a cyclic convolution.");
skip();
item("The operation count is $6K N \log N/2$.");
skip();
item("Implicit padding uses half the memory of explicit padding in higher dimensions as well.");

//title("Implicit Convolutions in Higher Dimensions");
//item("The transformed arrays are multiplied:");
//figure("cyrc_2dm","height=5cm");
//item("Once we have $F_kG_k$, we take the inverse transform to get $f*g$:");
//figure("cyrc_2dinv","height=4cm");

title("Alternatives");
item("The memory savings could be achieved more simply by using conventional padded transforms.");
item("This requires copying data, which is slow.");
step();
//skip();
item("Pruning: note that half of the FFTs in the $x$-direction are on zero-data.");
item("We can skip such transforms:");
figure("cyrc_prune","height=4cm");  // FIXME: indexedfigure
step();
remark("This is actually slower for large data sets due to memory-striding issues.");

title("Implicit Padding in Higher Dimensions");
item("Implicit padding is faster in two dimensions:");
figure("timing2c","height=15cm"); // TODO: redo (with N instead of m?)

title("Implicit Padding in Higher Dimensions");
item("The algorithm is easily extended to three dimensions:");
figure("timing3c","height=15cm"); // TODO: redo (with N instead of m?)

title("Centered Hermitian Data");
item("The Fourier transform of $f$ is \emph{centered} if ");
equation("\mathcal{F}(f)=\{F_k\}_{k=-N/2}^{N/2}.");
item("If $\{f_n\}_{n=0}^{N-1}$ is real-valued, then $\mathcal{F}(f)$ is \emph{Hermitian}:");
equation("F_{-k}=\conj{F}_k");
//item("Real-to-complex FFTs take $K\frac{N}{2}\log\frac{N}{2}$ multiplies.");
item("Padding centered data increases the array length by $50\%$");
figure("cyrc_23","height=4cm");
item("Phase-shifting is slower than explicit padding for centered data.");
item("Implicit padding uses 2/3 the memory of explicit padding.");

title("Hermitian Data");
item("The 1D implicit convolution is comparable to explicit padding:");
figure("timing1r","height=15cm"); // TODO: redo (with N instead of m?)

title("Hermitian Data");
item("And faster in higher dimensions:");
figure("timing2r","height=15cm"); // TODO: redo (with N instead of m?)


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
item("They use less memory and are faster than explicit zero-padding or phase-shift dealiasing."); // FIXME: check this!
item("Expanding into discontiguous arrays makes for easier programming.");
item("``Efficient Dealiased Convolutions without Padding\" submitted to SIAM Journal on Scientific Computing.");
item("A {\tt C++} implementation under the LGPL is available at {\tt http://fftwpp.sourceforge.net/}");
item("Uses SIMD routines when compiled with the Intel compiler.");
item("Uses the Fastest Fourier Transform in the West ({\tt http://fftw.org/}) for sub-transforms.");
figure("fftw-logo-med.eps","height=1.4cm"); 



if(false){
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
