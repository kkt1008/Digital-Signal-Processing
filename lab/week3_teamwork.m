%a
n=[0:39]; w = linspace(-pi,pi,1000);

h_n = -(0.6).^abs(n);

H = h_n * exp(-j*n.*w);

plot(w,H);


%b h(n)=sinc(0.2n)[u(n)-u(n-40)]
n=[0:39]; 
w = linspace(-pi,pi,1000);
h_n = sinc(0.2*n);

H_w = zeros(size(w));
for k=1:length(w)
    H_w(k)=sum(h_n .* exp(-j*w(k) *n ));
end

% magnitude와 phase
magnitude = abs(H_w);
phase = angle(H_w);

subplot(3,1,1);
plot(w,magnitude);
title('magnitude')

subplot(3,1,2);
plot(w,phase);
title('phase');

% 그래프
subplot(3,1,3);
plot(w,abs(H_w));
xlabel('omega'); ylabel('H');
title('H(e^jw) ');
