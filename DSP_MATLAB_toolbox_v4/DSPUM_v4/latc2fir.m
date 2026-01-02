function [a] = latc2fir(K,a0)
% Lattice form to FIR Direct form Conversion
% ------------------------------------------
% [a] = latc2fir(K,b0)
%  a = FIR direct form coefficients (prediction coefficients)
%  K = Lattice filter coefficients (reflection coefficients)
% a0 = Overall gain if \= 1 (optional)
%
if nargin == 1
    a0 = 1;
end
p = length(K);
B = 1; A = 1;
for m=1:1:p
	A = [A,0]+conv([0,K(m)],B);
	B = fliplr(A);
end
a = a0*A;
end