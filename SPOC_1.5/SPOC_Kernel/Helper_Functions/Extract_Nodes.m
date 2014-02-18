%=========================================================================
function [time_nodes, parameter_values, pdf_values] = ...
                                    Extract_Nodes(Discretization, Methods)
                                
                                
global CONSTANTS OFFLINE_TRAJECTORIES ...
       MESHED_PDF_WEIGHTS...
       DISCRETIZATION_VALUES MESHED_DISCRETIZATION_VALUES ...
       DIFFERENTIATION_MATRICES ...
       INTEGRATION_WEIGHTS MESHED_INTEGRATION_WEIGHTS

% Must run problem or set paths and run Problem_Definitions function before using.  
% After paths and problems definitions are set, function can be called for
% any method or discretization.

Calculate_Methods(Discretization, Methods);
Build_PDF();

time_nodes = DISCRETIZATION_VALUES{1};
parameter_values = MESHED_DISCRETIZATION_VALUES; 
pdf_values = MESHED_PDF_WEIGHTS;