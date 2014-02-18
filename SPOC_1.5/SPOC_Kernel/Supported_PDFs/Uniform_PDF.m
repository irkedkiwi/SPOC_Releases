% © 2012, CLAIRE WALTON. Some Rights Reserved.
%=========================================================================
function PDF = Uniform_PDF(discretization_values,parameter_index)

global CONSTANTS


w0 = CONSTANTS.ParameterSpace.W0(parameter_index);
wf = CONSTANTS.ParameterSpace.WF(parameter_index);

PDF = ones(length(discretization_values),1)./abs(wf-w0);