% © 2012, CLAIRE WALTON. Some Rights Reserved.

function Build_PDF()

global CONSTANTS OFFLINE_TRAJECTORIES ...
       MESHED_PDF_WEIGHTS...
       DISCRETIZATION_VALUES MESHED_DISCRETIZATION_VALUES ...
       DIFFERENTIATION_MATRICES ...
       INTEGRATION_WEIGHTS MESHED_INTEGRATION_WEIGHTS


MESHED_PDF_WEIGHTS = ones(length(MESHED_INTEGRATION_WEIGHTS),1);

if strcmp(CONSTANTS.ParameterSpace.PFD_Choices(1),'Independent')
    for i = 1:CONSTANTS.ParameterSpace.Dimension
        MESHED_PDF_WEIGHTS = ...
            MESHED_PDF_WEIGHTS.*feval(str2func(CONSTANTS.ParameterSpace.PFD_Choices{i+1}),MESHED_DISCRETIZATION_VALUES(:,i),i);
    end 
else
    MESHED_PDF_WEIGHTS = feval(str2func(CONSTANTS.ParameterSpace.PFD_Choices{2}),MESHED_DISCRETIZATION_VALUES);

    end
end
