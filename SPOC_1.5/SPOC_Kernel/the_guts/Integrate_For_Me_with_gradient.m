function [f,g] = Integrate_For_Me_with_gradient(inner_function, outer_function, inner_function_gradient, outer_function_gradient, controls, states)

global MESHED_DISCRETIZATION_VALUES ...
       MESHED_PDF_WEIGHTS...
       INTEGRATION_WEIGHTS...
       MESHED_INTEGRATION_WEIGHTS...
       DIFFERENTIATION_MATRICES...
       CONSTANTS
   
 
N = CONSTANTS.N;

g_u = zeros(N,CONSTANTS.Nu);
g_x = zeros(N,CONSTANTS.Nx);
f = 0;
for mesh_index = 1:size(MESHED_DISCRETIZATION_VALUES,1) 
    z  = inner_function(controls, states, mesh_index);
    
    [runningcost_gradient_controls, runningcost_gradient_states] =...
                inner_function_gradient(controls, states, mesh_index);
    
    if strcmpi(CONSTANTS.End_Cost, 'yes')==1
        [endcost_gradient_controls, endcost_gradient_states] =... 
            feval(str2func(CONSTANTS.End_Cost_Gradient), controls(end, :), states(end,:), mesh_index);
    else
        endcost_gradient_controls = zeros(1, CONSTANTS.Nu);
        endcost_gradient_states = zeros(1, CONSTANTS.Nx);
    end
    
    gradient_controls = runningcost_gradient_controls;
    gradient_states = runningcost_gradient_states;
    
    gradient_controls(end, :) = gradient_controls(end, :) + endcost_gradient_controls;
    gradient_states(end, :) = gradient_states(end, :) + endcost_gradient_states;

    time_integral = INTEGRATION_WEIGHTS{1}'*z;
    
    G_prime = outer_function_gradient(time_integral);
    
    weight = MESHED_INTEGRATION_WEIGHTS(mesh_index)*MESHED_PDF_WEIGHTS(mesh_index);
    
    u_partials = weight*G_prime*...
      repmat(INTEGRATION_WEIGHTS{1},1,CONSTANTS.Nu).*gradient_controls;
    x_partials = weight*G_prime*...
      repmat(INTEGRATION_WEIGHTS{1},1,CONSTANTS.Nx).*gradient_states;
    g_u = g_u + u_partials;
    g_x = g_x + x_partials;
    if strcmpi(CONSTANTS.End_Cost, 'no')==1
        f = f + weight*outer_function(time_integral) ; 
    else
        f = f + weight*(outer_function(time_integral)+feval(str2func(CONSTANTS.End_Cost), controls(end, :), states(end, :), mesh_index));
    end
end

temp_g=[g_u(:);g_x(:)];
[~, col] = find(CONSTANTS.indicators_first_row);
g = temp_g(col);

