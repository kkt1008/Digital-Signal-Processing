% 1조 문제

% 2조 문제
b1=[1,4,6,4,1]; a1=[1]; % b:분자, a:분모
b2=[1,0.5]; a2=[1,-1.7,1.36,-0.576];
b3=[1,0,1]; a3=[1,-0.2,1.53,-0.648];
% direct-->cascade c:계수, B:분자, A:분모
[c1,B1,A1]=dir2cas(b1,a1);
[c2,B2,A2]=dir2cas(b2,a2);
[c3,B3,A3]=dir2cas(b3,a3);

% direct-->parallel
[c1,B1,A1]=dir2par(b1,a1);
[c2,B2,A2]=dir2par(b2,a2);
[c3,B3,A3]=dir2par(b3,a3);

% 3조 문제
clear all; clc;
b=[1,-2.818,3.97,-2.818,1]; 
a=[1,-2.536,3.215,-2.054,0.656];

[c1,B1,A1]=dir2cas(b,a);
[c2,B2,A2]=dir2par(b,a);

% 4조 문제
% iir-->direct을 위해 conv함수를 이용하여 통분 진행
clear all; clc;
b1=[2,0,2]; a1=[1,-0.8,0.64];
b2=[2,-1,0]; a2=[1,-0.75,0];
b3=[1,2,1]; a3=[1,0,0.81];
% 분모 통분 먼저 진행
A=conv(conv(a1,a2),a3); % 최종 분모 통분값
% 분자 통분 항 별로 진행
B=conv(conv(a1,a2),b1)+conv(conv(a1,a3),b2)+conv(conv(a1,a2),b3);

% 8조 문제
clear all; clc;
b1=[0.5,2,1.5]; a1=[1,1,0.92];
b2=[1,3]; a2=[1,-1,0.8];
b3=[1,2,1]; a3=[1,0.5,0.5,-0.4];
b4=[3,-0.5,2]; a4=[1,0.4,0.4];

B1=conv(b1,b2);
A1=conv(a1,a2);
B2=conv(b3,b4);
A2=conv(a3,a4);

B_final=conv(A2,B1)+conv(A1,B2);
A_final=conv(A1,A2);

% 9조 문제---여러항을 하나의 배열표현으로 묶을 수 있음
clear all; clc;
c1=2; b1=[0.2,-0.3;0.4,0.5]; a1=[1,0.9,0.9;1,-0.8,0.8];
c2=1; b2=[2,1,-1; 3,4,5]; a2=[1,1.7,0.72;1,-1.5,0.56];
[B1,A1]=par2dir(c1,b1,a1);
[B2,A2]=cas2dir(c2,b2,a2);

B_final=conv(A2,B1)+conv(A1,B2);
A_final=conv(A1,A2);

% 10조 문제---frequency sampling form 변환
clear all; clc;
b=[exp((-0.9)*3),exp((-0.9)*2),exp((-0.9)*1),exp((-0.9)*0),exp((-0.9)*1),exp((-0.9)*2),exp((-0.9)*3)];
[C,B,A]=dir2fs(b);


% week6--teamwork 11조문제 
b=[2,3,5,-3,0,4,0,8 ,-7,4];
[b0,B,A]=dir2cas(b,1);
disp(b0); disp(B); disp(A);

h=[2,3,5,-3,0,4,0,8 ,-7,4];
[C,B,A]=dir2fs(h);
disp(C); disp(B); disp(A);

