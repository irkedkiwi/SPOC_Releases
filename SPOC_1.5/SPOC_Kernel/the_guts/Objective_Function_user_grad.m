function [F, G] = Objective_Function_user_grad(x)

global CONSTANTS...
       DIFFERENTIATION_MATRICES ...
       MESHED_DISCRETIZATION_VALUES ...
       MESHED_PDF_WEIGHTS...
       INTEGRATION_WEIGHTS...
       MESHED_INTEGRATION_WEIGHTS

 
DN = DIFFERENTIATION_MATRICES{1};
N_active_time_points = CONSTANTS.N_active_time_points;
Nx = CONSTANTS.Nx;
Nu = CONSTANTS.Nu;
N = CONSTANTS.N;

[controls, states] = Extract_States(x); 

[f,g] = Integrate_For_Me_with_gradient(str2func(CONSTANTS.Inner_Function), ...
                                       str2func(CONSTANTS.Outer_Function), ...
                                       str2func(CONSTANTS.Inner_Function_Gradient),...
                                       str2func(CONSTANTS.Outer_Function_Gradient), controls, states);
                                   
[nonlinear_dynamics, gradient_controls,gradient_states] ...
                                      = feval(str2func(CONSTANTS.Nonlinear_Dynamics),controls, states);
                              
[controls_indicator_vectors, states_indicator_vectors] ...
                                      = feval(str2func(CONSTANTS.Nonlinear_Dynamics_Sparsity));  
                                                                                 
dynamics_constraints = zeros(N_active_time_points,Nx);


for k=1:Nx
    if states_indicator_vectors{k}(k) == 1
        dynamics_constraints(:,k) = DN*states(:,k) - nonlinear_dynamics(1:N_active_time_points,k);
     
    else
        dynamics_constraints(:,k) = - nonlinear_dynamics(1:N_active_time_points,k);
     
    end
end

F = [f; dynamics_constraints(:)];

for k=1:Nx
    gradient_controls{k} = -1*gradient_controls{k}(1:N_active_time_points,:);
    gradient_states{k} = -1*gradient_states{k}(1:N_active_time_points,:);
end

constraint_rows=zeros(N_active_time_points*Nx, Nu*N+N*Nx);

for k=1:Nx
    
    down_chunk = N_active_time_points;
    over_chunk = N;
    
    control_indices = find(controls_indicator_vectors{k});
    for i = 1:length(control_indices)
        constraint_rows(down_chunk*(k-1)+1:down_chunk*k,...
                   over_chunk*(control_indices(i)-1)+1:over_chunk*(control_indices(i)-1)+N_active_time_points)...
            = diag(gradient_controls{k}(:, control_indices(i)));   
    end
    
    state_indices = find(states_indicator_vectors{k});
    for i = 1:length(state_indices)
        constraint_rows(down_chunk*(k-1)+1:down_chunk*k,...
                   over_chunk*(Nu+state_indices(i)-1)+1:over_chunk*(Nu+state_indices(i)-1)+N_active_time_points)...
            = diag(gradient_states{k}(:, state_indices(i)));   
    end

    if states_indicator_vectors{k}(k) == 1
        constraint_rows(down_chunk*(k-1)+1:down_chunk*k,...
                        over_chunk*(Nu+k-1)+1:over_chunk*(Nu+k-1)+N)...
        = DN + constraint_rows(down_chunk*(k-1)+1:down_chunk*k,...
                        over_chunk*(Nu+k-1)+1:over_chunk*(Nu+k-1)+N);     
    end       
end


big_G = [CONSTANTS.indicators_first_row; constraint_rows];

G_spots = find(CONSTANTS.indicators_first_row);

big_G(1,G_spots) = g(:); 
 
G =  big_G(CONSTANTS.G_indices);
 

