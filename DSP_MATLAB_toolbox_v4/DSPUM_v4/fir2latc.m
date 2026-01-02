function [K,a0] = fir2latc(a)
% FIR Direct form to All-Zero Lattice form Conversion
% ---------------------------------------------------
% [K,a0] = fir2latc(b)
%   K = Lattice filter coefficients (reflection coefficients)
%  a0 = First coefficient (or gain) of A(z), useful if \= 1
%   a = FIR direct form coefficients (prediction coefficients)
p = length(a)-1;
K = zeros(1,p);
a0 = a(1);
if a0 == 0
	error('a(1) is equal to zero')
end
A = a/a0;
for m=p+1:-1:2
	K(m-1) = A(m);
	B = fliplr(A);
	A = (A-K(m-1)*B)/(1-K(m-1)*K(m-1));
	A = A(1:m-1);
end
end