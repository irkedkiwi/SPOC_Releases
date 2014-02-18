% © 2012, CLAIRE WALTON. Some Rights Reserved.
%=========================================================================

function PDF = Beta_PDF(discretization_values, parameter_index)
%Returns column vector of pdf values at discretization_values

global CONSTANTS

w0 = CONSTANTS.ParameterSpace.W0(parameter_index);
wf = CONSTANTS.ParameterSpace.WF(parameter_index);
alpha = CONSTANTS.PDF.alpha(parameter_index);  
beta = CONSTANTS.PDF.beta(parameter_index);

% have to tranform range from [w0,wf] to [0,1] for matlab's betapdf
% Thus point is evaluated at (x-w0)/(wf-wf)

PDF = zeros(length(discretization_values),1);

for i = 1:length(discretization_values)
    x=discretization_values(i);
    PDF(i,1)=betapdf((x-w0)/(wf-w0),alpha,beta)/(wf-w0);
end