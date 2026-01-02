% p1-a
b1=[1,-1,-4,4]; a1=[1,-11/4,13/8,-1/4];
[R1,p1,C1]=residuez(b1,a1);

%p1-b
b2=[2,-7/4,-31/8,4]; a2=[1,-11/4,13/8,-1/4];
[R2,p2,C2]=residuez(b2,a2);

%p2
[delta, n]=impseq(0,0,10);

X1_check = -10*(0.5).^n+27*(0.25).^n-16*delta;
X1=filter(b1,a1,delta);
disp(X1_check-X1);

subplot(2,2,1);
stem(X1); title('X1');

subplot(2,2,2);
stem(X1_check); title('X1 check');

X2_check = 2.^n -10*0.5.^n + 27*(0.25).^n-16*delta;
X2=filter(b2,a2,delta);
disp(X2_check-X2);

subplot(2,2,3);
stem(X2); title('X2');

subplot(2,2,4);
stem(X2_check); title('X2 check');



%% P3
[H1,w1] = freqz(b1,a1,100);
magH1 = abs(H1); phaH1 = angle(H1);
[H2,w2] = freqz(b2,a2,100);
magH2 = abs(H2); phaH2 = angle(H2);

figure()
subplot(2,2,1);
plot(w1,magH1);
xlabel('Frequency [rad/sample]');
ylabel('|X_1(e^{j\omega})|');
grid on;
title('Magnitude Response of (a)');

subplot(2,2,2);
plot(w1,phaH1);
xlabel('Frequency [rad/sample]');
ylabel('Phase [rad]');
grid on;
title('Phase Response of (a)');

subplot(2,2,3);
plot(w2,magH2);
xlabel('Frequency [rad/sample]');
ylabel('|X_2(e^{j\omega})|');
grid on;
title('Magnitude Response of (b)');

subplot(2,2,4);
plot(w2,phaH2);
xlabel('Frequency [rad/sample]');
ylabel('Phase [rad]');
grid on;
title('Phase Response of (b)');
