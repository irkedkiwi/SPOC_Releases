% © 2012, CLAIRE WALTON. Some Rights Reserved.
%=========================================================================


function [x,F,xmul,Fmul,INFO, run_time] = Run_Optimization(Problem_Info)

global CONSTANTS



[Objective_lower, Objective_upper, epsilon, Max_u, Min_u, ...
          Max_x,Min_x,x_0] = feval(str2func(Problem_Info.Optimization_Bounds));

[Dynamics_lower, Dynamics_upper] = AutoCreate_Dynamics_Bounds(epsilon);
[Variables_lower, Variables_upper] = AutoCreate_Variable_Bounds(Min_u, Max_u, Min_x, Max_x, x_0);

xlow = Variables_lower;
xupp = Variables_upper;
Flow = [Objective_lower; Dynamics_lower];
Fupp = [Objective_upper; Dynamics_upper];
nF = length(Flow);
n = length(xlow);
[uguess, xguess] = feval(str2func(Problem_Info.Initial_Guess));
x = Package_States(uguess, xguess);
xstate = zeros(n,1);
ObjAdd = 0; 
ObjRow = 1;
xmul   =     zeros(n,1);
Fstate =   zeros(nF,1);
Fmul   =   zeros(nF,1);
Start = 1; 

CONSTANTS.Inner_Function = Problem_Info.Inner_Function;   
CONSTANTS.Outer_Function = Problem_Info.Outer_Function;
    
if strcmpi(CONSTANTS.End_Cost, 'yes')==1
    CONSTANTS.End_Cost = Problem_Info.End_Cost;
    if strcmpi(CONSTANTS.User_Gradient, 'yes')==1
        CONSTANTS.End_Cost_Gradient = Problem_Info.End_Cost_Gradient;
        CONSTANTS.End_Cost_Sparsity = Problem_Info.End_Cost_Sparsity;
    end
end

if strcmpi(CONSTANTS.User_Gradient, 'no')==1
    
    CONSTANTS.State_Dynamics = Problem_Info.State_Dynamics;
    
    Objective_Function = 'Objective_Function_auto_grad';
    disp('--------------------------------------------------------');
    disp('  ... ... ... ... ... ... ... ... ... ... ... ... ...   ');
    disp(' ... ... ... ... ... running snJac ... ... ... ... ...  ');
    disp('  ... ... ... ... ... ... ... ... ... ... ... ... ...   '); 
    
    [A,iAfun,jAvar,iGfun,jGvar] = snJac(Objective_Function,x,xlow,xupp,nF); 
    snseti ('Derivative option', 0);
    
elseif strcmpi(CONSTANTS.User_Gradient,'yes')==1 
    
    CONSTANTS.Inner_Function_Gradient = Problem_Info.Inner_Function_Gradient;
    CONSTANTS.Inner_Function_Sparsity = Problem_Info.Inner_Function_Sparsity;
    CONSTANTS.Outer_Function_Gradient = Problem_Info.Outer_Function_Gradient;
    CONSTANTS.Nonlinear_Dynamics = Problem_Info.Nonlinear_Dynamics;
    CONSTANTS.Nonlinear_Dynamics_Sparsity = Problem_Info.Nonlinear_Dynamics_Sparsity;
    CONSTANTS.Linear_Dynamics = Problem_Info.Linear_Dynamics;
    
    Objective_Function = 'Objective_Function_user_grad';
    [iGfun,jGvar] = S_Nonlinear_Gradient_Sparsity();
    %();
    [iAfun,jAvar,A] = S_Linear_Gradient();
    snseti ('Derivative option', 1);
    
else
   error('Error: Ivalid User_Gradient option!');
end
   

Optimization_Settings;

run_time = cputime;

disp('--------------------------------------------------------');
disp('  ... ... ... ... ... ... ... ... ... ... ... ... ...   ');
disp(' ... ... ... ... ... running snopt ... ... ... ... ...  ');
disp('  ... ... ... ... ... ... ... ... ... ... ... ... ...   '); 

[x,F,xmul,Fmul,INFO] = snoptcmex(Start, x, xlow, xupp, xmul, xstate, Flow,...
                                 Fupp, Fmul, Fstate, ObjAdd, ObjRow, A,...
                                 iAfun, jAvar, iGfun, jGvar,...
                                 Objective_Function);

snprintfile('off')                                   
run_time = cputime - run_time;
disp('--------------------------------------------------------');
disp(' Run Time:');
disp(run_time);
disp('--------------------------------------------------------');

