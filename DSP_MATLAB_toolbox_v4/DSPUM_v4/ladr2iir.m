function [b,a] = ladr2iir(K,beta)
% Lattice/Ladder form to IIR Direct form Conversion
% -------------------------------------------------
% [b,a] = ladr2iir(K,beta)
%     b = numerator polynomial coefficients
%     a = denominator polymonial coefficients
%     K = Lattice coefficients (reflection coefficients)
%  beta = Ladder coefficients
%
p = length(K); q = length(beta);
beta = [beta, zeros(1,p-q+1)];
J = 1; a = 1; A = zeros(p,p);
for m=1:1:p
    a = [a,0]+conv([0,K(m)],J);
    A(m,1:m) = -a(2:m+1);  J = fliplr(a);
end
b(p+1) = beta(p+1);
for m = p:-1:1
    A(m,1:m) = A(m,1:m)*beta(m+1);
    b(m) = beta(m) - sum(diag(A(m:p,1:p-m+1)));
end