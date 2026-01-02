% Kaiser 창과 Hann 창을 사용한 선형 위상 FIR 노치 필터 설계
clear all; close all;

% pass Band 영역
wp1 = 0;
wp2 = 0.25 * pi;

wp3 = 0.75 * pi;
wp4 = pi;
% stop Band 영역
ws1 = 0.35 * pi;
ws2 = 0.65 * pi;

% transition Band 영역
tr_width1 = ws1 - wp2; % 0.35pi - 0.25pi
tr_width2 = wp3 - ws2; % 0.75pi - 0.65pi
tr_width = min(tr_width1,tr_width2);

delta_1 = 0.025; % passband tolerance
delta_2 = 0.005; % stopband ripple

% stop band attenuation---필터길이 N값 결정
As = -20*log10(delta_2/(1+delta_1)); % As ~=4 6.2351
% pass band ripple
Rp = -20*log10((1-delta_1)/(1+delta_1));

% 필터 차수 
M0 = (As-7.95)/(14.36*tr_width/(2*pi));
%M = ceil(M0+1)+1; 
M= 61
% time step
n=0:1:M-1;

% cut-off frequency
wc1 = (0.25*pi + 0.35*pi)/2; % passband상한과 stopband하한의 평균
wc2 = (0.65*pi + 0.75*pi)/2; % passband하한과 stopband상한의 평균

% kaiser filter beta (As ~= 46)
beta = 0.5842*(As-21).^0.4+0.07886*(As-21);

% ideal BPF
hd1 = ideal_lp(pi,M); % wc=pi
hd2 = ideal_lp(wc2,M); % wc=wc2
hd3 = ideal_lp(wc1,M); % wc=w3
hd = hd1 - hd2 + hd3;

% kaiser filter
w_kai = kaiser(M,beta);

% kaiser 윈도우 적용
h = hd.*w_kai';
[db1,mag,pha,grd,w]=freqz_m(h,1);

delta_w = 2*pi/1000;

Rp1 = -(min(db1(wp2/delta_w+1:1:pi/delta_w)));
As1 = -round(max(db1(ws1/delta_w+1:1:ws2/delta_w)));

figure(1);
subplot(2,2,1); stem(n,hd); title('hd(n)');
subplot(2,2,2); stem(n,w_kai); title('w(n) kaiser window (M=56)');
subplot(2,2,3); stem(n,h); title('h(n)');
subplot(2,2,4); plot(w/pi,db1); title('Hr(w)'); axis([0 1 -100 10]);

% p2---hamming filter

M2 = (6.2*pi/tr_width) + 1;
n2 = 0:1:M2-1;
w_hann = hann(M2);
hd2 = ideal_lp(wc1,M2) + ideal_lp(pi,M2) -ideal_lp(wc2,M2);
h2 = hd2.*w_hann';
[db2,mag2,pha2,gra2,w2] = freqz_m(h2,1);
Rp2 = -min(db2(0/delta_w+1:1:wp1/delta_w));
As2 = -round(max(db2(ws1/delta_w+1:1:ws2/delta_w)));

figure(2);
subplot(2,2,1); stem(n2,hd2); title('hd(n)');
subplot(2,2,2); stem(n2,w_hann); title('w(n) hamming window');
subplot(2,2,3); stem(n2,h2); title('h(n)');
subplot(2,2,4); plot(w2/pi,db2); title('Hr(w)'); axis([0 1 -100 10]);
