clear; close all;

% Filter Design
I = 2; D = 5; Rp = 0.1; As = 30;
wxp = 0.2*pi; 
wxs1 = 0.45*pi; wxs2 = 0.5*pi;
% Filter parameters
wp = wxp/I; 
ws = min((2*pi/I-wxs1/I),(2*pi/D-wxs2/I));


[delta1,delta2] = db2delta(Rp,As);
f = [wp,ws]/pi;
[N,f,A,weights] = ...
  firpmord(f,[1,0],[delta1,delta2],2);
h = firpm(N,f,A,weights);

n = 0:length(h)-1;
[Hr,w,a,~] = hr_type1(h);
Hr_min = min(Hr);
H = abs(freqz(h,1,w));
Hdb = 20*log10(H/max(H));
min_attn = Hdb(Hr==Hr_min);

% Print results
fprintf('Stopband atten. = %4.2f (dB)\n',min_attn);

[Hr,omega,P,~] = ampl_res(h);
M = length(h);


% Plots
figure(1); 
subplot(3,1,1);
Hs = stem(0:1:M-1,h,'filled','MarkerSize',3);
grid on;
title('Impulse Response');
xlim([0,M-1]);
xlabel('n'); 
ylabel('h(n)');
set(gca,'XTickMode','manual','XTick',(0:10:M-1));

subplot(3,1,2);
Hdb2 = Hdb(w>ws);
ind = islocalmax(Hdb2);
min_attn = max(Hdb2(ind));
plot(w/pi,Hdb,'linewidth',1);
grid on;
title('Log-Magnitude Response');
axis([0,1,-70,10]); 
xlabel('frequency in \pi units'); 
ylabel('Magnitude (dB)');
set(gca,'XTickMode','manual','XTick',f);
set(gca,'YTickMode','manual','YTick',[-60,round(min_attn),0]);

subplot(3,1,3);
plot(omega/pi,Hr,'linewidth',1);
grid on;
title('Amplitude Response');
axis([0 1 -0.2 1.2]); 
xlabel('frequency in \pi units'); 
ylabel('Amplitude');
set(gca,'XTickMode','manual','XTick',f);
set(gca,'YTickMode','manual','YTick',[0,1]);

%% 3번
n = 0:500;   
x = sin(0.35*pi*n) + 2*cos(0.45*pi*n);
%y = upfirdn(x,h,1,D);

[delta1,delta2] = db2delta(Rp,As);
f = [wp,ws]/pi;
[N,f,A,weights] = ...
  firpmord(f,[1,0],[delta1,delta2],2);
h = firpm(N,f,A,weights);

y=upfirdn(x,h,1,D);