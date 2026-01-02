function [Px,xc] = pdf1(x,xmin,xmax,M)
% pdf1: Normalized Histogram as 1-D Probability Density Function (pdf)
%   [Px,x] = pdf1(xn,xmin,xmax,M)
%       Px: Normalized histogram over the range [xmin, xmax]
%        x: Histogram bib centers
%       xn: Data (observation) values
%     xmin: Minimun range value
%     xmax: Maximum range value
%        M: Number of bins

N = length(x);                         % observation count
edges = linspace(xmin,xmax,M);         % histogram boundaries
Dx = (xmax-xmin)/(M-1);                % Delta_x
xc = [xmin,edges(1:M-1)+Dx/2,xmax];    % Bin centers
edges = [-inf,edges,inf];              % Augment boundaries
Nx = histc(x,edges); Nx = Nx(1:end-1); % Histogram
Px = Nx/(N*Dx);                        % Normalized Histogram
end