function [K,beta] = iir2ladr(b,a)
% IIR Direct form to pole-zero Lattice/Ladder form Conversion
% -----------------------------------------------------------
% [K,beta] = iir2ladr(b,a)
%     K = Lattice coefficients (reflection coefficients), [K1,...,KN]
%  beta = Ladder Coefficients, [C0,...,CN]
%     b = Numerator polynomial coefficients (deg <= Num deg)
%     a = Denominator polynomial coefficients
%
a1 = a(1); a = a/a1; b = b/a1;
q = length(b); p = length(a);
if q > p
     error('   *** length of b must be <= length of a ***')
end
b = [b, zeros(1,p-q)]; K = zeros(1,p-1);
A = zeros(p-1,p-1); beta = b;
for m = p-1:-1:1
     A(m,1:m) = -a(2:m+1)*beta(m+1);
     K(m) = a(m+1);  J = fliplr(a);
     a = (a-K(m)*J)/(1-K(m)*K(m));  a = a(1:m);
     beta(m) = b(m) + sum(diag(A(m:p-1,1:p-m)));
end