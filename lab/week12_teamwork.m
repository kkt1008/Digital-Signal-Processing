%% p1
clc; clear all;
% Digital Filter Specifications:
wp = 0.8*pi;                         % digital Passband freq in Hz
ws = 0.95*pi;                         % digital Stopband freq in Hz
Rp = 0.5;                              % Passband ripple in dB
As = 45;                             % Stopband attenuation in dB

% Analog Prototype Specifications: Inverse mapping for frequencies
T = 1; Fs = 1/T;                     % Set T=1
OmegaP = (2/T)*tan(wp/2);            % Prewarp Prototype Passband freq
OmegaS = (2/T)*tan(ws/2);            % Prewarp Prototype Stopband freq
ep = sqrt(10^(Rp/10)-1);             % Passband Ripple parameter
Ripple = sqrt(1/(1+ep*ep));          % Passband Ripple
Attn = 1/(10^(As/20));               % Stopband Attenuation

% Analog Chebyshev Prototype Filter Calculation:
[cs,ds] = afd_chb1(OmegaP,OmegaS,Rp,As);
%%*** Chebyshev-1 Filter Order =  4 

% Bilinear transformation:
[b,a] = bilinear(cs,ds,T);
[C,B,A] = dir2cas(b,a)

% Plotting
figure;
[db,mag,pha,grd,w] = freqz_m(b,a);
subplot(2,2,1); plot(w/pi,mag,'linewidth',1); axis([0,1,0,1.1]); 
title('Magnitude Response','verticalalignment','baseline');
xlabel('frequency in \pi units','verticalalignment','middle'); 
ylabel('|H|','verticalalignment','baseline'); %fontchar;
set(gca,'XTickMode','manual','XTick',[0,0.2,0.3,1]);
set(gca,'YTickmode','manual','YTick',[0,Attn,Ripple,1]); grid;

subplot(2,2,3); plot(w/pi,db,'linewidth',1); axis([0,1,-40,5]);
title('Magnitude in dB','verticalalignment','baseline');
xlabel('frequency in \pi units','verticalalignment','middle'); 
ylabel('Decibels','verticalalignment','baseline'); %fontchar;
set(gca,'XTickMode','manual','XTick',[0,0.2,0.3,1]);
set(gca,'YTickmode','manual','YTick',[-50,-15,-1,0]); grid;
set(gca,'YTickLabelMode','manual','YTickLabel',['50';'15';' 1';' 0']);

subplot(2,2,2); plot(w/pi,pha/pi,'linewidth',1); axis([0,1,-1,1]);
title('Phase Response','verticalalignment','baseline');
xlabel('frequency in \pi units','verticalalignment','middle'); 
ylabel('Radians in \pi units','verticalalignment','baseline'); %fontchar;
set(gca,'XTickMode','manual','XTick',[0,0.2,0.3,1]);
set(gca,'YTickmode','manual','YTick',[-1,0,1]); grid;

n=0:30;
h=impz(b,a,n);
subplot(2,2,4); stem(n,h,'linewidth',1); axis([0,30,-2,2]);
title('Impulse Response in Digital Filter','verticalalignment','baseline');
xlabel('time in seconds','verticalalignment','middle'); 
ylabel('Samples','verticalalignment','baseline'); %fontchar;

figure;
[ha,x,t]=impulse(b,a); % 아날로그 필터의 임펄스 응답 계산
subplot(1,2,1); stem(n,h,'linewidth',1);
title('Digital filter-Impulse Response');

subplot(1,2,2); plot(t,ha,'linewidth',1);
title('Analog filter-Impulse Response');

%% p2
clc; clear all;

% Digital Filter Specifications:
wp = 0.8*pi;                         % digital Passband freq in Hz
ws = 0.95*pi;                         % digital Stopband freq in Hz
Rp = 0.5;                              % Passband ripple in dB
As = 45;                             % Stopband attenuation in dB

% Analog Prototype Specifications: Inverse mapping for frequencies
T = 1; Fs = 1/T;                     % Set T=1
OmegaP = (2/T)*tan(wp/2);            % Prewarp Prototype Passband freq
OmegaS = (2/T)*tan(ws/2);            % Prewarp Prototype Stopband freq
ep = sqrt(10^(Rp/10)-1);             % Passband Ripple parameter
Ripple = sqrt(1/(1+ep*ep));          % Passband Ripple
Attn = 1/(10^(As/20));               % Stopband Attenuation

% Analog Butterworth Prototype Filter Calculation:
[cs,ds] = afd_butt(OmegaP,OmegaS,Rp,As);
%%*** Butterworth Filter Order =  6 

% Bilinear transformation:
[b,a] = bilinear(cs,ds,T);
[C,B,A] = dir2cas(b,a)

%% Plotting
figure;
[db,mag,pha,grd,w] = freqz_m(b,a); axis([0,1,0,1.1]);
subplot(2,2,1); plot(w/pi,mag,'linewidth',1); 
title('Magnitude Response','verticalalignment','baseline');
xlabel('frequency in \pi units','verticalalignment','middle'); 
ylabel('|H|','verticalalignment','baseline'); %fontchar;
set(gca,'XTickMode','manual','XTick',[0,0.2,0.3,1]);
set(gca,'YTickmode','manual','YTick',[0,Attn,Ripple,1]); grid;

subplot(2,2,3); plot(w/pi,db,'linewidth',1); axis([0,1,-40,5]);
title('Magnitude in dB','verticalalignment','baseline');
xlabel('frequency in \pi units','verticalalignment','middle'); 
ylabel('Decibels','verticalalignment','baseline'); %fontchar;
set(gca,'XTickMode','manual','XTick',[0,0.2,0.3,1]);
set(gca,'YTickmode','manual','YTick',[-50,-15,-1,0]); grid;
set(gca,'YTickLabelMode','manual','YTickLabel',['50';'15';' 1';' 0']);

subplot(2,2,2); plot(w/pi,pha/pi,'linewidth',1); axis([0,1,-1,1]);
title('Phase Response','verticalalignment','baseline');
xlabel('frequency in \pi units','verticalalignment','middle'); 
ylabel('Radians in \pi units','verticalalignment','baseline'); %fontchar;
set(gca,'XTickMode','manual','XTick',[0,0.2,0.3,1]);
set(gca,'YTickmode','manual','YTick',[-1,0,1]); grid;

n=0:30;
h=impz(b,a,n);
subplot(2,2,4); stem(n,h,'linewidth',1); axis([0,30,-2,2]);
title('Impulse Response in Digital filter','verticalalignment','baseline');
xlabel('time in seconds','verticalalignment','middle'); 
ylabel('Samples','verticalalignment','baseline'); %fontchar;
set(gca,'XTickMode','manual','XTick',[0,0.2,0.3,1]);
set(gca,'YTickmode','manual','YTick',[0:2:10]); grid;

%% p3
clc; clear all;

% Digital Filter Specifications:
wp = 0.8*pi;                         % digital Passband freq in Hz
ws = 0.95*pi;                         % digital Stopband freq in Hz
Rp = 0.5;                              % Passband ripple in dB
As = 45;                             % Stopband attenuation in dB

% Analog Prototype Specifications: Inverse mapping for frequencies
T = 1; Fs = 1/T;                     % Set T=1
OmegaP = (2/T)*tan(wp/2);            % Prewarp Prototype Passband freq
OmegaS = (2/T)*tan(ws/2);            % Prewarp Prototype Stopband freq
ep = sqrt(10^(Rp/10)-1);             % Passband Ripple parameter
Ripple = sqrt(1/(1+ep*ep));          % Passband Ripple
Attn = 1/(10^(As/20));               % Stopband Attenuation

% Analog Chebyshev-2 Prototype Filter Calculation:
[cs,ds] = afd_chb2(OmegaP,OmegaS,Rp,As);
%%*** Chebyshev-2 Filter Order =  4 

% Bilinear transformation:
[b,a] = bilinear(cs,ds,T);
[C,B,A] = dir2cas(b,a)

% Plotting
figure;
[db,mag,pha,grd,w] = freqz_m(b,a);
subplot(2,2,1); plot(w/pi,mag,'linewidth',1); axis([0,1,0,1.1]); 
title('Magnitude Response','verticalalignment','baseline');
xlabel('frequency in \pi units','verticalalignment','middle'); 
ylabel('|H|','verticalalignment','baseline'); %fontchar;
set(gca,'XTickMode','manual','XTick',[0,0.2,0.3,1]);
set(gca,'YTickmode','manual','YTick',[0,Attn,Ripple,1]); grid;

subplot(2,2,3); plot(w/pi,db,'linewidth',1); axis([0,1,-40,5]); 
title('Magnitude in dB','verticalalignment','baseline');
xlabel('frequency in \pi units','verticalalignment','middle'); 
ylabel('Decibels','verticalalignment','baseline'); %fontchar;
set(gca,'XTickMode','manual','XTick',[0,0.2,0.3,1]);
set(gca,'YTickmode','manual','YTick',[-50,-15,-1,0]); grid;
set(gca,'YTickLabelMode','manual','YTickLabel',['50';'15';' 1';' 0']);

subplot(2,2,2); plot(w/pi,pha/pi,'linewidth',1); axis([0,1,-1,1]); 
title('Phase Response','verticalalignment','baseline');
xlabel('frequency in \pi units','verticalalignment','middle'); 
ylabel('Radians in \pi units','verticalalignment','baseline'); %fontchar;
set(gca,'XTickMode','manual','XTick',[0,0.2,0.3,1]);
set(gca,'YTickmode','manual','YTick',[-1,0,1]); grid;

n=0:30;
h=impz(b,a,n);
subplot(2,2,4); stem(n,h,'linewidth',1); axis([0,30,-2,2]);
title('Impulse Response in Digital filter','verticalalignment','baseline');
xlabel('Time in seconds','verticalalignment','middle'); 
ylabel('Samples','verticalalignment','baseline'); %fontchar;
set(gca,'XTickMode','manual','XTick',[0,0.2,0.3,1]);
set(gca,'YTickmode','manual','YTick',[0:2:10]); grid;

%% p4
clc; clear all;
% Digital Filter Specifications:
wp = 0.8*pi;                         % digital Passband freq in Hz
ws = 0.95*pi;                         % digital Stopband freq in Hz
Rp = 0.5;                              % Passband ripple in dB
As = 45;                             % Stopband attenuation in dB

% Analog Prototype Specifications: Inverse mapping for frequencies
T = 1; Fs = 1/T;                     % Set T=1
OmegaP = (2/T)*tan(wp/2);            % Prewarp Prototype Passband freq
OmegaS = (2/T)*tan(ws/2);            % Prewarp Prototype Stopband freq
ep = sqrt(10^(Rp/10)-1);             % Passband Ripple parameter
Ripple = sqrt(1/(1+ep*ep));          % Passband Ripple
Attn = 1/(10^(As/20));               % Stopband Attenuation

% Analog Elliptic Prototype Filter Calculation:
[cs,ds] = afd_elip(OmegaP,OmegaS,Rp,As);
%%*** Elliptic Filter Order =  3 

% Bilinear transformation:
[b,a] = bilinear(cs,ds,T);
[C,B,A] = dir2cas(b,a)


% Plotting
figure;
[db,mag,pha,grd,w] = freqz_m(b,a);
subplot(2,2,1); plot(w/pi,mag,'linewidth',1); axis([0,1,0,1.1]);
title('Magnitude Response','verticalalignment','baseline');
xlabel('frequency in \pi units','verticalalignment','middle');
ylabel('|H|','verticalalignment','baseline'); %fontchar;
set(gca,'XTickMode','manual','XTick',[0,0.2,0.3,1]);
set(gca,'YTickmode','manual','YTick',[0,Attn,Ripple,1]); grid;

subplot(2,2,3); plot(w/pi,db,'linewidth',1); axis([0,1,-40,5]);
title('Magnitude in dB','verticalalignment','baseline');
xlabel('frequency in \pi units','verticalalignment','middle');
ylabel('Decibels','verticalalignment','baseline'); %fontchar;
set(gca,'XTickMode','manual','XTick',[0,0.2,0.3,1]);
set(gca,'YTickmode','manual','YTick',[-50,-15,-1,0]); grid;
set(gca,'YTickLabelMode','manual','YTickLabel',['50';'15';' 1';' 0']);

subplot(2,2,2); plot(w/pi,pha/pi,'linewidth',1); axis([0,1,-1,1]);
title('Phase Response','verticalalignment','baseline');
xlabel('frequency in \pi units','verticalalignment','middle');
ylabel('Radians in \pi units','verticalalignment','baseline'); %fontchar;
set(gca,'XTickMode','manual','XTick',[0,0.2,0.3,1]);
set(gca,'YTickmode','manual','YTick',[-1,0,1]); grid;

n=0:30;
h=impz(b,a,n);
subplot(2,2,4); stem(n,h,'linewidth',1); axis([0,30,-2,2]);
title('Impulse Response in Digital filter','verticalalignment','baseline');
xlabel('frequency in \pi units','verticalalignment','middle');
ylabel('Samples','verticalalignment','baseline'); %fontchar;
set(gca,'XTickMode','manual','XTick',[0,0.2,0.3,1]);
set(gca,'YTickmode','manual','YTick',[0:5:15]); grid;

