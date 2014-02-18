function [controls, states] = Extract_States(x)


global CONSTANTS DIFFERENTIATION_MATRICES

N = size(DIFFERENTIATION_MATRICES{1},2);
Nx =  CONSTANTS.Nx;

Nu = CONSTANTS.Nu;

    
controls = zeros(N,Nu);
states = zeros(N,Nx);
controls(:)=x(1:N*Nu);
states(:)=x(N*Nu+1:end);
