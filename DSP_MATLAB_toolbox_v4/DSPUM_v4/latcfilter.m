function [y] = latcfilter(K,x,a0)
% LATTICE form realization of FIR filters
% ---------------------------------------
% y = latcfilter(K,x,a0)
%   y = output sequence
%   K = LATTICE filter (reflection) coefficient array
%   x = input sequence
%  a0 = Overall Gain (optional)
%
if nargin == 1
    a0 = 1;
end
Nx = length(x)-1;
x = a0*x;
p = length(K); %K = K(2:M+1);
fg = [x; [0 x(1:Nx)]];
for m = 1:p
    fg = [1,K(m);K(m),1]*fg;
    fg(2,:) = [0 fg(2,1:Nx)];
end
y = fg(1,:);
end