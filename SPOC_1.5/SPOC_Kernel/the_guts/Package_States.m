function [x] = Package_States(controls, states)

global CONSTANTS DIFFERENTIATION_MATRICES


N = size(DIFFERENTIATION_MATRICES{1},2);
  
Nx =  CONSTANTS.Nx;

Nu = CONSTANTS.Nu;

x = [controls(:);states(:)];