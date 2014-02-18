function [Variables_lower, Variables_upper] =...
        AutoCreate_Variable_Bounds(Min_u, Max_u, Min_x, Max_x, x_0)
    
global CONSTANTS
    
N = CONSTANTS.N;
Nu = CONSTANTS.Nu;
Nx = CONSTANTS.Nx;

u_lower = repmat(Min_u,N,1);
u_upper = repmat(Max_u,N,1);

x_lower = repmat(Min_x,N,1);
x_upper = repmat(Max_x,N,1);

for k=1:Nx
    x_lower(1,k) = x_0(k);
    x_upper(1,k) = x_0(k);
end

Variables_lower = Package_States(u_lower, x_lower);
Variables_upper = Package_States(u_upper, x_upper);