function [F] = Objective_Function_auto_grad(x)

global CONSTANTS...
       DIFFERENTIATION_MATRICES ...
       MESHED_DISCRETIZATION_VALUES ...
       MESHED_PDF_WEIGHTS...
       INTEGRATION_WEIGHTS...
       MESHED_INTEGRATION_WEIGHTS
DN = DIFFERENTIATION_MATRICES{1};
N_active_time_points = size(DIFFERENTIATION_MATRICES{1},1);

Nx = CONSTANTS.Nx;
    

 
[controls, states] = Extract_States(x); 

f = Integrate_For_Me(str2func(CONSTANTS.Inner_Function),str2func(CONSTANTS.Outer_Function), controls, states);

dynamics = feval(str2func(CONSTANTS.State_Dynamics),controls, states);
   
dynamics_constraints = zeros(N_active_time_points,Nx); 
for k=1:Nx
    dynamics_constraints(:,k) = DN*states(:,k) - dynamics(1:N_active_time_points,k);
end

constraint_vector = Package_Constraints(dynamics_constraints);

    
F = [f;constraint_vector];
