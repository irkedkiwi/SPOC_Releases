function [iAfun,jAvar,A] = S_Linear_Gradient()

global CONSTANTS DIFFERENTIATION_MATRICES

[gradient_controls, gradient_states,...
 controls_indicator_vectors, states_indicator_vectors] ...
                    = feval(str2func(CONSTANTS.Linear_Dynamics()));
[nonlinear_controls_indicator_vectors, nonlinear_states_indicator_vectors] ...
                    = feval(str2func(CONSTANTS.Nonlinear_Dynamics_Sparsity()));                                

DN = DIFFERENTIATION_MATRICES{1};

N = CONSTANTS.N;
N_active_time_points = CONSTANTS.N_active_time_points;

Nx = CONSTANTS.Nx;
Nu = CONSTANTS.Nu;
    

for k=1:Nx
    gradient_controls{k} = -1*gradient_controls{k}(1:N_active_time_points,:);
    gradient_states{k} = -1*gradient_states{k}(1:N_active_time_points,:);
end

obj_row = zeros(1,Nu*N+Nx*N);
constraint_rows=zeros(N_active_time_points*Nx, Nu*N+N*Nx);

for k=1:Nx
    
    down_chunk = N_active_time_points;
    over_chunk = N;
    
    control_indices=find(controls_indicator_vectors{k});
    for i = 1:length(control_indices)
        constraint_rows(down_chunk*(k-1)+1:down_chunk*k,...
                       over_chunk*(control_indices(i)-1)+1:over_chunk*(control_indices(i)-1)+N_active_time_points)...
                    = diag(gradient_controls{k}(:, control_indices(i)));   
    end
    
    state_indices=find(states_indicator_vectors{k});
    
    for i = 1:length(state_indices)
        constraint_rows(down_chunk*(k-1)+1:down_chunk*k,...
                        over_chunk*(Nu+state_indices(i)-1)+1:over_chunk*(Nu+state_indices(i)-1)+N_active_time_points)...
            = diag(gradient_states{k}(:, state_indices(i)));   
    end

    if nonlinear_states_indicator_vectors{k}(k) == 0
        constraint_rows(down_chunk*(k-1)+1:down_chunk*k,...
                        over_chunk*(Nu+k-1)+1:over_chunk*(Nu+k-1)+N)...
        = DN + constraint_rows(down_chunk*(k-1)+1:down_chunk*k,...
                        over_chunk*(Nu+k-1)+1:over_chunk*(Nu+k-1)+N);     
    end
end


[iAfun,jAvar,A] = find(A);




