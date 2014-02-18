% © 2012, CLAIRE WALTON. Some Rights Reserved.
%=========================================================================
% FUNCTION [Dn, x, w] = euldiff(N)
% =================================
% INPUT = Number of node points, N
%
% OUTPUT =  Euler Diff Matrix, Dn
%           Equispaced points, x
%           Integration Weights, w
% This function calculates the differentiation matrix Dn that is
% obtained by using forward Euler's method on n equally spaced points.
% These values are calculated over the interval [-1,1] and need to be
% transformed to the interval of the problem
% Note that this does not return a square matrix. It returns an (N-1)xN
% matrix
%=========================================================================
function [Dn,x,w]=euldiff(N)

h=2/(N-1);
x = (-1:h:1)';
w = ones(N,1).*h;
w(N,1)=0;
v1=ones(1,N-1);
v2=-1*ones(1, N);
temp=(diag(v1,1)+diag(v2))./h;
Dn = temp(1:N-1,:);