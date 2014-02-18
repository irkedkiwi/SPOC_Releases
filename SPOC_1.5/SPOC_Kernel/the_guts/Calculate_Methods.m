% © 2012, CLAIRE WALTON. Some Rights Reserved.

function Calculate_Methods(Discretization, Methods)

global CONSTANTS OFFLINE_TRAJECTORIES ...
       MESHED_PDF_WEIGHTS...
       DISCRETIZATION_VALUES MESHED_DISCRETIZATION_VALUES ...
       DIFFERENTIATION_MATRICES ...
       INTEGRATION_WEIGHTS MESHED_INTEGRATION_WEIGHTS



if Methods(1)==0
    disp('--------------------------------------------------------')
    disp('Methods:');
    disp('Time: Euler');
    [Dn, Time_Array, Weights]=euldiff(Discretization(1));
end

if Methods(1)==1
    disp('--------------------------------------------------------')
    disp('Methods:');
    disp('Time: Pseudospectral with Legendre Points');
    [Dn, Time_Array, Weights]=LGdiff(Discretization(1));
      
end
if Methods(1)>1
    disp('--------------------------------------------------------')
    disp('Method does not exist for time domain!')
    return;
end

T0=CONSTANTS.Time.T0;
TF=CONSTANTS.Time.TF;
DISCRETIZATION_VALUES{1} = .5*(TF+T0)+.5*(TF-T0).*Time_Array;

DIFFERENTIATION_MATRICES{1} = (2/(TF-T0)).*Dn;

INTEGRATION_WEIGHTS{1} = .5*(TF-T0).*Weights;


Size=CONSTANTS.ParameterSpace.Dimension+1;

if Methods(2)<=1 
    for i=2:Size

        if Methods(i)==0
            str=['Parameter ',num2str(i-1), ': Euler'];
            disp(str);
            [Dn, Parameter_Array, Weights]=euldiff(Discretization(i));
        end   

        if Methods(i)==1
            str=['Parameter ',num2str(i-1), ': Pseudospectral with Legendre Points'];
            disp(str);
            [Dn, Parameter_Array, Weights]=LGdiff(Discretization(i));
            % Get differentiation matrix, legendre points, and quadrature weights
        end
        if Methods(i)==2
        disp('--------------------------------------------------------')
        disp('Cannot call sparse grid for partial parameter space!')
        return;
        end
        W0=CONSTANTS.ParameterSpace.W0(i-1);
        WF=CONSTANTS.ParameterSpace.WF(i-1);
        DISCRETIZATION_VALUES{i} = .5*(WF+W0)+.5*(WF-W0).*Parameter_Array;

        DIFFERENTIATION_MATRICES{i} = (2/(WF-W0)).*Dn;

        INTEGRATION_WEIGHTS{i} = .5*(WF-W0).*Weights;
        disp('--------------------------------------------------------') 
    end

    MESHED_DISCRETIZATION_VALUES = DISCRETIZATION_VALUES{2}(:);
    MESHED_INTEGRATION_WEIGHTS = INTEGRATION_WEIGHTS{2}(:);
    for i=3:Size
        chunk_length = size(MESHED_DISCRETIZATION_VALUES,1);
        temp1 = repmat(MESHED_DISCRETIZATION_VALUES,length(DISCRETIZATION_VALUES{i}),1);
        temp2 = zeros(size(temp1,1),1);
        index = 1;
        for j = 1:length(DISCRETIZATION_VALUES{i})
            for k = 1:chunk_length
                temp2(index) = DISCRETIZATION_VALUES{i}(j);
                index = index+1;
            end
        end
        MESHED_DISCRETIZATION_VALUES = [temp1, temp2];

        temp1 = repmat(MESHED_INTEGRATION_WEIGHTS,length(INTEGRATION_WEIGHTS{i}),1);
        temp2 = zeros(size(temp1,1),1);
        index = 1;
        for j = 1:length(INTEGRATION_WEIGHTS{i})
            for k = 1:chunk_length
                temp2(index) = INTEGRATION_WEIGHTS{i}(j);
                index = index+1;
            end
        end
        MESHED_INTEGRATION_WEIGHTS = temp1.*temp2;
    end

else
    if Methods(2)>2
    disp('--------------------------------------------------------')
    disp('Method does not exist!')
    return;
    
    elseif Methods(2)==2   
        str='Parameter Space: Sparse grid with univariate nested quadrature rules as basis - delayed Kronrod-Patterson rules, see Knut Petras (2003)';
        disp(str);
        
        [grid, Weights]=nwspgr('KPU', CONSTANTS.ParameterSpace.Dimension, Discretization(2));
  
        MESHED_DISCRETIZATION_VALUES = grid;
        MESHED_INTEGRATION_WEIGHTS = Weights;

        for i=2:Size
            W0=CONSTANTS.ParameterSpace.W0(i-1);
            WF=CONSTANTS.ParameterSpace.WF(i-1);
            MESHED_DISCRETIZATION_VALUES(:,i-1) = W0+(WF-W0).*MESHED_DISCRETIZATION_VALUES(:,i-1);
            MESHED_INTEGRATION_WEIGHTS = abs(WF-W0).*MESHED_INTEGRATION_WEIGHTS; %this needs to be checked  
        end
    end
end
    


