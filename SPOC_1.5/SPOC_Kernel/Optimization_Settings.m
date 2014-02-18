%=========================================================================
% For modifying snopt settings.
%
% Almost any setting can be modified here, except gradient settings and 
% settings which modify the loaction of the cost function relative to other 
% variables.
%=========================================================================
global CONSTANTS
% Example settings: 
temporary_superbasics_limit = 1.25*CONSTANTS.N*(CONSTANTS.Nx+CONSTANTS.Nu);
    % For the general case, I would advise not changing this. If you're
    % getting INFO=33., increase that multiplier by using the line of code below.
snseti('Superbasics limit', temporary_superbasics_limit);

snseti ('Major print level',10^3);
snseti ('Minor print level',10^3);
snseti ('Scale option',1); %default is 1. 0 is recommended if values are smallish, say less than 100
snseti ('Verify level',3);
snseti ('Major Iteration limit', 10^5);
snseti ('Minor Iteration limit', 10^4);
snseti ('Iterations limit', 10^5);
snsetr ('Major feasibility tolerance', 10^(-4))
snsetr ('Minor feasibility tolerance', 10^(-3))
snsetr ('Major optimality tolerance', 10^(-4))
snsetr ('Minor optimality tolerance', 10^(-3))
 
 fdirout = './Results';
 filename = [fdirout, '/', 'snopt_printout', datestr(clock), '.txt'];

 snscreen('off')
 snprintfile( filename )
