function [iGfun,jGvar] = S_Nonlinear_Gradient_Sparsity()

global CONSTANTS DIFFERENTIATION_MATRICES

[controls_indicator_vectors, states_indicator_vectors] ...
                = feval(str2func(CONSTANTS.Nonlinear_Dynamics_Sparsity));                        

DN = DIFFERENTIATION_MATRICES{1};
N = CONSTANTS.N;
N_active_time_points = CONSTANTS.N_active_time_points;

Nx = CONSTANTS.Nx;
Nu = CONSTANTS.Nu;

    
indicator_G_u = cell(1,CONSTANTS.Nx);
indicator_G_x = cell(1,CONSTANTS.Nx);
for k=1:Nx
    indicator_G_u{k} = repmat(controls_indicator_vectors{k},N_active_time_points,1);
    indicator_G_x{k} = repmat(states_indicator_vectors{k},N_active_time_points,1);
end    


[obj_indicator_vector_controls, obj_indicator_vector_states] = feval(str2func(CONSTANTS.Inner_Function_Sparsity));

if strcmpi(CONSTANTS.End_Cost, 'yes')==1
    [endcost_indicator_vector_controls, endcost_indicator_vector_states] = feval(str2func(CONSTANTS.End_Cost_Sparsity));
else
    endcost_indicator_vector_controls = zeros(1, CONSTANTS.Nu);
    endcost_indicator_vector_states = zeros(1, CONSTANTS.Nx);
end

temp_obj = [repmat(obj_indicator_vector_controls,N,1),...
        repmat(obj_indicator_vector_states,N,1)];
    
temp_endcost = zeros(N, CONSTANTS.Nu + CONSTANTS.Nx);
temp_endcost(end, :) = [endcost_indicator_vector_controls, endcost_indicator_vector_states];

temp = temp_obj+temp_endcost; 
CONSTANTS.indicators_first_row = temp(:)';

CONSTANTS.N_first_row_elements = length(find(temp));

obj_row = CONSTANTS.indicators_first_row; 

constraint_rows=zeros(N_active_time_points*Nx, Nu*N+N*Nx);

for k=1:Nx
    
    down_chunk = N_active_time_points;
    over_chunk = N;
    
    control_indices=find(controls_indicator_vectors{k});
    for i = 1:length(control_indices)
        constraint_rows(down_chunk*(k-1)+1:down_chunk*k,...
                over_chunk*(control_indices(i)-1)+1:over_chunk*(control_indices(i)-1)+N_active_time_points)...
            = diag(indicator_G_u{k}(:, control_indices(i)));  

    end
    
    state_indices=find(states_indicator_vectors{k});
    
    for i = 1:length(state_indices)
        constraint_rows(down_chunk*(k-1)+1:down_chunk*k,...
                over_chunk*(Nu+state_indices(i)-1)+1:over_chunk*(Nu+state_indices(i)-1)+N_active_time_points)...
            = diag(indicator_G_x{k}(:, state_indices(i)));  
    end

    if states_indicator_vectors{k}(k) == 1
        constraint_rows(down_chunk*(k-1)+1:down_chunk*k,...
                        over_chunk*(Nu+k-1)+1:over_chunk*(Nu+k-1)+N)...
        = DN + constraint_rows(down_chunk*(k-1)+1:down_chunk*k,...
                        over_chunk*(Nu+k-1)+1:over_chunk*(Nu+k-1)+N);  

    end        
end
 
indicator_G = [obj_row; constraint_rows];
size(indicator_G) 

[iGfun,jGvar,G] = find(indicator_G);

CONSTANTS.iGfun = iGfun;
CONSTANTS.jGvar = jGvar;

CONSTANTS.G_indices = sub2ind(size(indicator_G),...
                        CONSTANTS.iGfun(:),CONSTANTS.jGvar(:));

    

