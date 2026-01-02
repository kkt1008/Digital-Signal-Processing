clear all; clc;
ws = 0.4*pi; wp = 0.6*pi;
As = 50; Rp = 0.004;
M = 37; alpha = (M-1)/2;
l = 0:M-1; wl = (2*pi/M)*l;
T1 = 1/3; T2 = 2/3;

% stopband & passband
wp1 = 0.6*pi; wp2 = 1*pi;
ws1 = 0; ws2 = 0.4*pi;

% P1: 주파수 샘플링 기법으로 HPF 설계
Hrs = [zeros(1,9), T1, T2, ones(1,16), T2, T1, zeros(1,8)];
Hdr = [0,0,1,1]; % ideal HPF
wdl = [0,0.4,0.6,1]; %T1,T2를 고려한 0~pi사이 정규화

k1 = 0:floor((M-1)/2);
k2 = floor((M-1)/2)+1:M-1;
angH = [-alpha*(2*pi)/M*k1, alpha*(2*pi)/M*(M-k2)];

H = Hrs.*exp(j*angH);
h = real(ifft(H,M));
[db,mag,pha,grd,w] = freqz_m(h,1);
[Hr,ww,a,L] = hr_type1(h);

% stopband & passband
wp1 = 0.6*pi; wp2 = 1*pi;
ws1 = 0; ws2 = 0.4*pi;

delta_ww = 2*pi/1000
ws_range = 1:1:ws/delta_ww+1

Rp = -min(db(wp1/delta_ww+1:1:wp2/delta_ww))
As = - round(max(db(ws_range)))

% 1. Impulse Response
figure(1)
stem(0:M-1, h, 'filled', 'LineWidth', 1.5);
title('Impulse Response h[n]');
xlabel('n'); ylabel('h[n]');
grid on;

% 2. Amplitude Response (from hr_type1)
figure(2);
plot(ww/pi, Hr, 'LineWidth', 1.5);
title('Amplitude Response H_r(\omega)');
xlabel('\omega / \pi'); ylabel('H_r(\omega)');
grid on;

% 3. Magnitude Response (in dB)
figure(3);
plot(w/pi, db, 'LineWidth', 1.5);
title('Magnitude Response (dB)');
xlabel('\omega / \pi'); ylabel('|H(\omega)| (dB)');
grid on;

% P2: T1, T2 변화시키며 A_s 계산
max_As = 0;

% max_As = -inf;  % As는 음수니까 초기값을 -inf로

for i = 0.05:0.05:1
    T2 = i;
    for j = 0.05:0.05:1
        T1 = j;

        % Hr 배열 새로 생성
        Hr_test = [zeros(1,9), T1, T2, ones(1,16), T2, T1, zeros(1,8)];
        
        % angH (기존 코드와 동일)
        H_test = Hr_test .* exp(1j * angH);
        h_test = real(ifft(H_test));
        
        % freqz 계산
        [db_test, mag_test, pha_test, grd_test, omega_test] = freqz_m(h_test, 1);
        
        % As 계산
        idx_ws2 = round(ws2 / delta_ww) + 1;
        ws_range = idx_ws2:1:501;
        As_meas = -max(db_test(ws_range));

        if As_meas > max_As
            max_As = As_meas;
            T1T2_opt = [j, i];
        end
    end
end

% 최종 결과 출력
fprintf('As = %.2f dB at T1 = %.2f, T2 = %.2f\n', max_As, T1T2_opt(1), T1T2_opt(2));
