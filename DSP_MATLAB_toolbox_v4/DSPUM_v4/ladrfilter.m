function [y] = ladrfilter(K,beta,x)
% LATTICE/LADDER form realization of IIR filters
% ----------------------------------------------
% [y] = ladrfilter(K,beta,x)
%    y = output sequence
%    K = LATTICE (reflection) coefficient array
% beta = LADDER coefficient array
%    x = input sequence
%
Nx = length(x); y = zeros(1,Nx);
p = length(beta); f = zeros(p,Nx); g = zeros(p,Nx+1);
f(p,:) = x;
for n = 2:1:Nx+1
   for m = p:-1:2
       f(m-1,n-1) = f(m,n-1) - K(m-1)*g(m-1,n-1);
       g(m,n) = K(m-1)*f(m-1,n-1) + g(m-1,n-1);
   end
   g(1,n) = f(1,n-1);
end
y = beta*g(:,2:Nx+1);