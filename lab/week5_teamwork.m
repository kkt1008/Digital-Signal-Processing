% 15조 문제
clear all; clc;
x1=[4,1,-1,4]; N1=4;
x2=[2,0,0,0,-1,0,0,0]; N2=8;
x3=[1,0,-1,0]; N3=4;
x4=[0,0,2j,0,2j,0]; N4=6;
% dfs함수를 이용해 Xk계수 구하기
X1=dfs(x1,N1);
X2=dfs(x2,N2);
X3=dfs(x3,N3);
X4=dfs(x4,N4);
% Xk계수로 주파수-위상 응답 구하기
magX1=abs(X1); angX1=angle(X1)*180/pi;
magX2=abs(X2); angX2=angle(X2)*180/pi;
magX3=abs(X3); angX3=angle(X3)*180/pi;
magX4=abs(X4); angX4=angle(X4)*180/pi;
% 주파수-위상 응답 그래프 결과
subplot(8,2,1); stem(magX1); 
subplot(8,2,2); stem(angX1);
subplot(8,2,3); stem(magX2);
subplot(8,2,4); stem(angX2);
subplot(8,2,5); stem(magX3);
subplot(8,2,6); stem(angX3);
subplot(8,2,7); stem(magX4);
subplot(8,2,8); stem(angX4);
% idfs함수를 이용해 Xk의 idfs계산-->다시 xn구함
x1=idfs(X1,N1);
x2=idfs(X2,N2);
x3=idfs(X3,N3);
x4=idfs(X4,N4);
%dftmtx 함수를 통해 변환행렬 Wn을 구할 수 있음
xn=[0,1,2,3]; Xk=dftmtx(4)*xn'; % 행렬곱셈을 위해 전치
Xn=dfs(xn,4);

% 13조 문제
% -80~80에서 각 주기 함수 그리기
clear all; clc;
N1=40; x1=zeros(1,N1);
x1(1:20)=5*sin(0.1*pi*(0:19)); % 특정 주기만 함수값 나타내는 표기법 주목!!
final_x1 = [x1,x1,x1,x1]; % -80~80사이에서 N=40인 주기가 4번 반복하므로 4번으로 설정

N2=80; x2=zeros(1,N2);
x2(1:20)=5*sin(0.1*pi*(0:19)); % 특정 주기만 함수값 나타내는 표기법 주목!!
final_x2 = [x2,x2];  % -80~80사이에서 N=80인 주기가 2번 반복하므로 2번으로 설정

N3=80; x3=zeros(1,N3);
x3(1:40)=5*sin(0.1*pi*(0:39)); % 특정 주기만 함수값 나타내는 표기법 주목!!
final_x3 = [x3,x3];  % -80~80사이에서 N=80인 주기가 2번 반복하므로 2번으로 설정

subplot(3,1,1); n1=[-80:79]; stem(n1,final_x1);
subplot(3,1,2); n2=[-80:79]; stem(n2,final_x2);
subplot(3,1,3);n3=[-80:79]; stem(n3,final_x3);

% 12조 문제
clear all; clc;
N=40; n=[0:N-1]; nunPeriods =4; 
x1 = zeros(1,N);
x1(1:20) = 5*sin(0.1*pi*(0:19));

x2 = zeros(1,N);
x2(1:20) = x1(1:20);
x2(21:40) = -x1(1:20);

x1_tilde = repmat(x1,1,nunPeriods); %배열을 반복해서 가로방향(1) 복사
x2_tilde = repmat(x2,1,nunPeriods);
subplot(2,1,1); stem(x1_tilde);
subplot(2,1,2); stem(x2_tilde);

% magnitude & phase
X1=dfs(x1(1:40),40);
X2=dfs(x2(1:40),40);
k = -40/2:40/2-1;
magX1=abs(X1); angX1=angle(X1);
magX2=abs(X2); angX2=angle(X2);

subplot(2,2,1); stem(k,magX1); title('magnitude of x1(k)');
subplot(2,2,2); stem(k,magX2); title('magnitude of x2(k)');
subplot(2,2,3); stem(k,angX1); title('phase of x1(k)');
subplot(2,2,4); stem(k,angX2); title('phase of x2(k)');

theoretical_X2 = zeros(1,N);
for kk =0:N-1
    theoretical_X2(kk+1) = X1(kk+1)*(1-exp(-1j*2*pi*20*kk/N));
end
error_X2 = X2-theoretical_X2;
error_X2_shifted = [error_X2(N/2+1:N), error_X2(1:N/2)];


% 11조 문제
clear all; clc;
N=100; x_n=[1:50 50:-1:1];
% DFS of x(n)
WN = dftmtx(N);
Xk = x_n*WN;

% frequency sampling으로 10-point IDFS를 이용하여 y(n)도출
N_sa2 =10;
w=2*pi*(0:N_sa2-1)/N_sa2;
diricfun = zeros(100,N_sa2);
for i=N_sa2
    diricfun(:,i)=diric(w(i)-2*pi*(0:N-1)/N,N).*exp(-1j*(w(i)-2*pi*(0:N-1)/N)*(N-1)/2);
end
% DTFT
x_dtft=Xk*diricfun;
% IDFS
WNcon = conj(dftmtx(N_sa2));
figure();
x_n2 = x_dtft*WNcon/N_sa2;
stem(0:N_sa2-1, x_n2);

% 10조 문제
clear all; clc;
n = [0:100]; N=101; k = [0:N-1]; w=2*pi*k/N;
x_n = sinc((n-50)/2).^2;
X_k = fft(x_n,N);

figure(1);
subplot(2,1,1); stem(k,abs(X_k)); title('magnitude');
subplot(2,1,2); stem(k,angle(X_k)); title('phase')

w=2*pi*(0:N-1)/N;
X_e = x_n*exp(-j*w.*n');
figure(2);
subplot(2,1,1); stem(w,abs(X_e)); title('magnitude');
subplot(2,1,2); stem(w,angle(X_e)); title('Phase');

figure(3);
subplot(2,1,1); stem(w,abs(X_k)); hold on;
stem(w,abs(X_e)); hold off;
subplot(2,1,2); stem(w,angle(X_k)); hold on;
stem(w,angle(X_e)); hold off;


% p1번
% 문제 P1: x(n) 계산 및 10-point DFT 계산
clear all; close all; clc;

% 주어진 시퀀스 X(k), k = 0, 1, ..., 5 (6-point DFT)
X = [10, -2+3j, 3+4j, 2-3j, 4+5j, 12, 4-5j, 2+3j, 3+4j,-2+3j];

% 10-point IDFT 계산 (X(k) -> x(n))
x = ifft(X, 10);
[x_even, x_odd]=circevod(x);

subplot(2,2,1);stem(x_even); title('x even part');
subplot(2,2,2); stem(x_odd); title('x odd part');

% 10-point DFT 계산
X = fft(x, 10);

[X_even, X_odd]=circevod(X);

subplot(2,2,3);stem(X_even); title('X even part');
subplot(2,2,4); stem(X_odd); title('X odd part');

% 결과 출력
disp('x(n) 값 (6-point IDFT 결과, 길이 6):');
disp(x_6point);

disp('x(n) 값 (10-point로 확장, 0 패딩):');
disp(x_10point);

disp('X(k) 값 (10-point DFT 결과):');
disp(X_10point);


% p2번
% 문제 P2: 변환된 시퀀스의 DFT 계산 (P1 결과 기반)
clear all; close all; clc;

% P1에서 계산된 x(n) (길이 10)
x_6point = [4.8333, 0.5-1.4434j, 0.1667+0.5774j, 0.8333, 0.1667-0.5774j, 0.5+1.4434j];
x = [x_6point, zeros(1, 4)]; % 길이 10으로 확장 (0 패딩)
N = 10;

% 1. x1(n) = x((2-n))_10
n = 0:N-1;
indices_x1 = mod(2 - n, N); % (2-n) mod 10
x1 = x(indices_x1 + 1); % MATLAB 인덱스는 1부터 시작
X1 = fft(x1, N);

% 2. x2(n) = x((n+5))_10
indices_x2 = mod(n + 5, N); % (n+5) mod 10
x2 = x(indices_x2 + 1);
X2 = fft(x2, N);

% 3. x3(n) = x(n) * x((-n))_10
indices_x3 = mod(-n, N); % (-n) mod 10
x3 = x .* x(indices_x3 + 1);
X3 = fft(x3, N);

% 4. x4(n) = x(n)_10 * x((-n))_10
% x(n)_10은 x(n) 자체 (이미 0 <= n <= 9)
x4 = x .* x(indices_x3 + 1); % 동일하게 계산됨 (x3과 동일)
X4 = fft(x4, N);

% 5. x5(n) = x(n) * e^(-j*4*pi*n/5)
x5 = x .* exp(-1j * 4 * pi * n / 5);
X5 = fft(x5, N);

% 결과 출력
disp('x(n) (길이 10):');
disp(x);

disp('x1(n) 및 X1(k):');
disp('x1(n) = '); disp(x1);
disp('X1(k) = '); disp(X1);

disp('x2(n) 및 X2(k):');
disp('x2(n) = '); disp(x2);
disp('X2(k) = '); disp(X2);

disp('x3(n) 및 X3(k):');
disp('x3(n) = '); disp(x3);
disp('X3(k) = '); disp(X3);

disp('x4(n) 및 X4(k):');
disp('x4(n) = '); disp(x4);
disp('X4(k) = '); disp(X4);

disp('x5(n) 및 X5(k):');
disp('x5(n) = '); disp(x5);
disp('X5(k) = '); disp(X5);

