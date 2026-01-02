% 2주차 조별 과제 
% 1-a번
n1=[-10:10];
x1=8*stepseq(0,-10,10)-5*stepseq(4,-10,10)-7*stepseq(6,-10,10)+6*stepseq(7,-10,10);

subplot(3,1,1);
stem(n1,x1); title('2주차 1-a번');
xlabel('n'); ylabel('x1(n)');

%1-b번
n2=[-10:10];
x2= exp(0.1.*n2).*(stepseq(-4,-10,10)-stepseq(8,-10,10));

[xe,xo,m]=evenodd(x2,n2);

subplot(3,1,2);
xlabel('m'); ylabel('x2(n)');
stem(m,xe); 

subplot(3,1,3);
xlabel('m'); ylabel('x2(n)');
stem(m,xo); 


%2번
% 2-1
n=[-10:15]; a1=2; a2=4;
new_x1 = sigshift(x1,n,2);
new_x2 = sigshift(x2,n,2);

y1=0.5*(a1*new_x1./0.5 + a2*new_x2./0.5);
y2=a1*new_x1 + a2*new_x2;

subplot(1,2,1);
title('2주차 2-a번');
xlabel('n'); ylabel('x(n)');
stem(n,y1); 

subplot(1,2,2);
xlabel('n'); ylabel('x2(n)');
stem(n,y2); 

% 2-2
