function m = mod1(n,N)
% Computes m = (n mod N) index
% ----------------------------
% m = mod1(n,N)
m = rem(n,N);
m = m+N;
m = rem(m,N);
