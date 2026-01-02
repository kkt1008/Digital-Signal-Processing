% DSP using MATLAB Toolbox (DSPUM_v4)
% Version 4 01-Oct-2015
%
% Files
%   afd_butt       - Analog Lowpass Filter Design: Butterworth
%   afd_chb1       - Analog Lowpass Filter Design: Chebyshev-1
%   afd_chb2       - Analog Lowpass Filter Design: Chebyshev-2
%   afd_elip       - Analog Lowpass Filter Design: Elliptic
%   cas2dir        - CASCADE-to-DIRECT form conversion
%   casfiltr       - CASCADE form realization of IIR and FIR filters
%   cheb1hpf       - IIR Highpass filter design using Chebyshev-1 prototype
%   circevod       - signal decomposition into circular-even and circular-odd parts
%   circonvt       - N-point circular convolution between x1 and x2: (time-domain)
%   cirshftt       - Circular shift of m samples wrt size N in sequence x: (time domain)
%   conv_m         - Modified convolution routine for signal processing
%   cplxcomp       - Compare two complex-valued pairs
%   dfs            - Computes Discrete Fourier Series Coefficients
%   dft            - Computes Discrete Fourier Transform
%   dir2cas        - DIRECT-form to CASCADE-form conversion (cplxpair version)
%   dir2fs         - Direct form to Frequency Sampling form conversion
%   dir2par        - DIRECT-form to PARALLEL-form conversion
%   evenodd        - Real signal decomposition into even and odd parts
%   fir2latc       - FIR Direct form to All-Zero Lattice form Conversion
%   freqs_m        - Computation of s-domain frequency response: Modified version
%   freqz_m        - Modified version of freqz subroutine
%   hr_type1       - Computes Amplitude response Hr(w) of a Type-1 LP FIR filter
%   hr_type2       - Computes Amplitude response of Type-2 LP FIR filter
%   hr_type3       - Computes Amplitude response Hr(w) of a Type-3 LP FIR filter
%   hr_type4       - Computes Amplitude response of Type-4 LP FIR filter
%   hsolpsav       - High-speed Overlap-Save method of block convolutions using FFT
%   ideal_lp       - Ideal LowPass filter computation
%   idfs           - Computes Inverse Discrete Fourier Series
%   idft           - Computes Inverse Discrete Transform
%   iir2ladr       - IIR Direct form to pole-zero Lattice/Ladder form Conversion
%   imp_invr       - Impulse Invariance Transformation from Analog to Digital Filter
%   impseq         - Generates x(n) = delta(n-n0); n1 <= n,n0 <= n2
%   ladr2iir       - Lattice/Ladder form to IIR Direct form Conversion
%   ladrfilter     - LATTICE/LADDER form realization of IIR filters
%   latc2fir       - Lattice form to FIR Direct form Conversion
%   latcfilter     - LATTICE form realization of FIR filters
%   lms            - Algorithm for Coefficient Adjustment
%   mod1           - Computes m = (n mod N) index
%   mulaw_c        - mu-law compressor
%   mulaw_e        - mu-law expander
%   OnesComplement - Sign-Magnitude format integer to b-bit One's-Complement format conversion
%   ovrlpsav       - Overlap-Save method of block convolution
%   par2dir        - PARALLEL-to-DIRECT form conversion
%   parfiltr       - PARALLEL form realization of IIR filters
%   pdf1           - pdf1: Normalized Histogram as 1-D Probability Density Function (pdf)
%   PSD            - Computation of PSD using Autocorrelation Lags
%   Q_Rounding     - Binary equivalent xq by rounding of x to B fraction bits
%   Q_Truncation   - Binary equivalent xq by truncation of x to B fraction bits
%   QCoeff         - Coefficient Quantization using N=L+B bit Representation with Rounding operation
%   QFix           - Fixed-point Arithmetic using (B+1)-bit Representation
%   randnMV        - randnMV: Multivariate Gaussian Random Vector Generator
%   sdir2cas       - DIRECT-form to CASCADE-form conversion in s-plane
%   sigadd         - implements y(n) = x1(n)+x2(n)
%   sigfold        - implements y(n) = x(-n)
%   sigmult        - implements y(n) = x1(n)*x2(n)
%   sigshift       - implements y(n) = x(n-n0)
%   StatModelR     - Statistical Model (Rounding) for A/D Quantization error and its Distribution
%   stepseq        - Generates x(n) = u(n-n0); n1 <= n,n0 <= n2
%   suptitle       - Puts a title above all subplots.
%   TwosComplement - Sign-Magnitude format integer to b-bit Two's-Complement format conversion
%   u_buttap       - Unnormalized Butterworth Analog Lowpass Filter Prototype
%   u_chb1ap       - Unnormalized Chebyshev-1 Analog Lowpass Filter Prototype
%   u_chb2ap       - Unnormalized Chebyshev-2 Analog Lowpass Filter Prototype
%   u_elipap       - Unnormalized Elliptic Analog Lowpass Filter Prototype
%   VarGain        - Computation of variance-gain for the output noise process
%   zmapping       - Frequency band Transformation from Z-domain to z-domain
