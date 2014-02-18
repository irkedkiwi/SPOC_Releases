function f = Integrate_For_Me(inner_function, outer_function, controls, states)



global CONSTANTS...
       DISCRETIZATION_VALUES...
       MESHED_DISCRETIZATION_VALUES ...
       MESHED_PDF_WEIGHTS...
       INTEGRATION_WEIGHTS...
       MESHED_INTEGRATION_WEIGHTS


f = 0;
for mesh_index = 1:size(MESHED_DISCRETIZATION_VALUES,1) %loop through all necessary permutations of parameter values 

    z  = inner_function(controls, states, mesh_index); 
    time_integral = INTEGRATION_WEIGHTS{1}'*z;
    
    weight = MESHED_INTEGRATION_WEIGHTS(mesh_index)*MESHED_PDF_WEIGHTS(mesh_index);
    
    if strcmpi(CONSTANTS.End_Cost, 'yes')==1
        f = f + weight*(outer_function(time_integral)+feval(str2func(CONSTANTS.End_Cost), controls(end, :), states(end, :), mesh_index));   
    else
        f = f + weight*outer_function(time_integral);
    end
end