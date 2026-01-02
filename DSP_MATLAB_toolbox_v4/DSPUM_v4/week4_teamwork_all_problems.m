% 1조 문제
% a
b=[1,0]; a=[1,-0.8];
[delta, n]=impseq(0,0,10);
X=filter(b,a,delta);
X=0.8^3*X;

stem(n,X); title('1번-a 결과');

X_check=(0.8)^3*(n>=3); 
disp(X-X_check);

% b
[delta,n]=impseq(0,0,10);
b1=[1,0]; a1=[1,0.5];
X1=filter(b1,a1,delta);

b2=[1,0]; a2=[1,-0.6];
X2=filter(b2,a2,delta);

X=X1+X2; stem(n,X); title('1번-b 결과')

% c
[delta,n]=impseq(0,0,10);
b1=[0,1]; a1=[2,1];
X1=filter(b1,a1,delta); 

b2=[0,1]; a2=[-1/0.6, 1];
X2=filter(b2,a2,delta); 

X=fliplr(X1+X2); stem(n,X);
title('1번-c 결과');

X_check= ((-0.5).^-n+ (0.6).^n).*(n<=10);
disp(X-X_check);

