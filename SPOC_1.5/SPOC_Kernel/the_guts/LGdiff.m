
function [Dn,x,w]=LGdiff(n);
% This function calculates the Legendre-Gauss nodes and weights.
% It also calculates the associated differentiation matrix from poldif
% from the differentiation suite.

ab=r_jacobi(n);
xw=gauss(n,ab);
% first column of xw has the nodes, and the 2nd column has the weights.
x=xw(:,1);
w=xw(:,2);
Dn=poldif(x,1);


