
function [Dynamics_lower, Dynamics_upper] = AutoCreate_Dynamics_Bounds(epsilon)

global CONSTANTS DIFFERENTIATION_MATRICES

N_active_time_points = size(DIFFERENTIATION_MATRICES{1},1);

N_dynamics_constraints = CONSTANTS.Nx;
    
Dynamics_lower = zeros(N_dynamics_constraints*N_active_time_points,1)-epsilon;
Dynamics_upper = zeros(N_dynamics_constraints*N_active_time_points,1)+epsilon; 
