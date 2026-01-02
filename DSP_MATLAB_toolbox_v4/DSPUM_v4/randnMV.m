function X = randnMV(N,mu,C)
% randnMV: Multivariate Gaussian Random Vector Generator
%   Generates N vectors in X given mean mu and covariance matrix C
%   mu should be a Mx1 column vector; C should be MxM
%   Generated X is NxM
% X = randnMV(N,mu,C)
% 
mu = mu(:); M = length(mu);
Y = randn(N,M);
X = sqrtm(C)*Y'+ mu*ones(1,N); X = X';
end

