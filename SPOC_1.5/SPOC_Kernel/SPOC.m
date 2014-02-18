% © 2013, CLAIRE WALTON. Some Rights Reserved.
% cwalton@soe.ucsc.edu
% irkedkiwi@gmail.com
%=========================================================================
%                               SPOC 1.4                                 % 
%=========================================================================
function [Results] = SPOC(Problem_Info, Discretization, Methods)

global CONSTANTS OFFLINE_TRAJECTORIES ...
       MESHED_PDF_WEIGHTS...
       DISCRETIZATION_VALUES MESHED_DISCRETIZATION_VALUES ...
       DIFFERENTIATION_MATRICES ...
       INTEGRATION_WEIGHTS MESHED_INTEGRATION_WEIGHTS


build_time = cputime;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% disp(Problem_Info.Name);
% disp('--------------------------------------------------------');
disp('Discretization:');
disp(Discretization);

%Set automatically:  
CONSTANTS.N = Discretization(1);



feval(str2func(Problem_Info.Problem_Definitions), Discretization, Methods);

Calculate_Methods(Discretization, Methods);
CONSTANTS.N_active_time_points = size(DIFFERENTIATION_MATRICES{1},1);

CONSTANTS.time_nodes = DISCRETIZATION_VALUES{1};    

Build_PDF();

build_time = cputime - build_time;
disp('--------------------------------------------------------');
disp(' Build Time:');
disp(build_time);
disp('--------------------------------------------------------');

[x,F,xmul,Fmul,INFO, run_time] = ...
                Run_Optimization(Problem_Info);

[Results] = Result_Outputs(x,F,xmul,Fmul);

Results.build_time = build_time;
Results.run_time = run_time;
Results.INFO = INFO;


