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
	  "John C. Bowman and Malcolm Roberts (University of Alberta)",
          date="April 13, 2010",
	  url="www.math.ualberta.ca/$\sim$bowman/talks");

/* Abstract:
                 The Fastest Convolution in the West
                         John C. Bowman
                      University of Alberta

Efficient algorithms have recently been developed for calculating dealiased
linear convolution sums without the expense of conventional zero-padding or
phase-shift techniques. For one-dimensional in-place convolutions, the
memory requirements are identical with the zero-padding technique, with the
important distinction that the additional work memory need not be
contiguous with the input data. This decoupling of data and work arrays
dramatically reduces the memory and computation time required to evaluate
higher-dimensional in-place convolutions. The memory savings is achieved by
computing the in-place Fourier transform of the data set in blocks, rather than
all at once. The technique also allows one to dealias the hyperconvolutions
that arise on Fourier transforming cubic and higher powers.

Implicitly dealiased convolutions can be built on top of state-of-the-art
adaptive fast Fourier transform libraries like FFTW. Vectorized
multidimensional implementations for the complex and centered Hermitian
(pseudospectral) cases have already been implemented in the open-source
software FFTW++. With the advent of this library, writing a high-performance
dealiased pseudospectral code for solving nonlinear partial differential
equations has now become a relatively straightforward exercise.

*/

newslide(stepping=false);
title("Outline",newslide=false);
item("Discrete Convolutions");
subitem("Cyclic vs.\ Linear");
subitem("Standard vs.\ Centered");
subitem("Complex vs.\ Hermitian");
item("Dealiasing");
subitem("Zero Padding");
subitem("Phase-shift dealiasing");
item("Implicit Padding in 1D, 2D, and 3D:");
subitem("Standard Complex");
subitem("Centered Hermitian");
subitem("Biconvolution");
item("Conclusions");

title("Discrete Convolutions");

item("Discrete linear convolution sums based on the fast Fourier transform
(FFT) algorithm~\cite{Gauss1866,Cooley65} have become important tools
for:");
subitem("image filtering;");
subitem("digital signal processing;");
subitem("correlation analysis;");
subitem("pseudospectral simulations.");

title("Discrete Cyclic Convolution");

item("The FFT provides an efficient tool for computing the {\it discrete cyclic convolution}");
equation("\sum_{p=0}^{N-1} F_p G_{k-p},");
remark("where the vectors $F$ and $G$ have period $N$.");

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

title("\mbox{Asymptote: 2D \& 3D Vector Graphics Language}");
//asyinclude("logo3");
center("Andy Hammerlindl, John C. Bowman, Tom Prince");
center("\tt http://asymptote.sf.net");
center("(freely available under the Lesser GNU Public License)");

title("Asymptote Lifts \TeX\ to 3D");
//asyinclude("label3solid");
skip();
center("\tt http://asymptote.sf.net");
center("Acknowledgements: Orest Shardt (U. Alberta)");


bibliography("refs");

//  LocalWords:  hypercitebracket dotover cdot mathop nolimits Im vl wl sim
//  LocalWords:  citename Dealiasing dealiasing asyinclude Hammerlindl Orest
//  LocalWords:  Shardt nocite
