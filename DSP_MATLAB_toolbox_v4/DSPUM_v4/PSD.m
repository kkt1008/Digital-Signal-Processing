function [Sx,omega] = PSD(Rx,maxlag,Nfft)
%PSD Computation of PSD using Autocorrelation Lags
%   [Sx,omega] = PSD(Rx,lags,Nfft)
% Sx: Computed PSD values
% omega: Digital frequency array in pi units from -1 to 1
% Rx: Autocorrelations from -maxlag to +maxlag
% maxlag: Maximum lag index (must be >= 10)
% Nfft: FFT size (must be >= 512)

Nfft2 = Nfft/2;
M = 2*maxlag+1; % Bartlett window length
Rx = bartlett(M).*Rx(:); % Windowed autocorrelations
Rxzp = [zeros(Nfft2-maxlag,1);Rx;zeros(Nfft2-maxlag-1,1)];
Rxzp = ifftshift(Rxzp); %Zero-padding and circular shifting
Sx = fftshift(real(fft(Rxzp))); % PSD
Sx = [Sx;Sx(1)]; % Circular shifting
omega = linspace(-1,1,Nfft+1); % Frquencis in units of pi
end