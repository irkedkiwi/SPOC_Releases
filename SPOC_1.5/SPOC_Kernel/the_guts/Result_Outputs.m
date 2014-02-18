% © 2012, CLAIRE WALTON. Some Rights Reserved.
function [Results] =...
            Result_Outputs(x,F,xmul,Fmul)
     
global CONSTANTS ...
       MESHED_PDF_WEIGHTS...
       DISCRETIZATION_VALUES...
       MESHED_DISCRETIZATION_VALUES ...
       DIFFERENTIATION_MATRICES ...
       INTEGRATION_WEIGHTS ...
       MESHED_INTEGRATION_WEIGHTS

[Results.controls, Results.states] = Extract_States(x);
Results.objective_value = F(1);
Results.time_nodes = DISCRETIZATION_VALUES{1};
Results.F = F;

    Results.constants = CONSTANTS; 
    Results.meshed_pdf_weights = MESHED_PDF_WEIGHTS;
    Results.discretization_values = DISCRETIZATION_VALUES; 
    Results.meshed_discretization_values = MESHED_DISCRETIZATION_VALUES;
    Results.differentiation_matrices = DIFFERENTIATION_MATRICES;
    Results.integration_weights = INTEGRATION_WEIGHTS;
    Results.meshed_integration_weights = MESHED_INTEGRATION_WEIGHTS;
