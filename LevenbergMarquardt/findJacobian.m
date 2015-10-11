%% Chantzi Efthymia - Optimization - Assignment 1 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates an approximation of the Jacobian using       %
% finite differences. It is called when the user has not supplied       %
% the gradient/Jacobian inside the objective function                   %
%                                                                       %
% %%%% Inputs %%%%                                                      %
% func: name of user's objective function                               %
% F: residual vector at x                                               %
% x: vector of unknown parameters                                       %
% t: vector of t data                                                   %
% y: vector of y data                                                   %
%                                                                       %
% %%%% Outputs %%%%                                                     %
% J: matrix of Jacobian (dimensions: rows_F x rows_x)                   %                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J = findJacobian(func, F, x, t, y)

numOfParam = length(x);                     %number of unknown parameters
residualRows = length(F);                   %number of rows of residual-vector function
J = zeros(residualRows, numOfParam);        %dimensions of the Jacobian

for i = 1 : numOfParam
   
    dx = 0.25*(10^-8);                 %very small change
    x_dif = x;
    x_dif(i) = x_dif(i) + dx;          %change in x vector with respect to the i parameter
    F_dif = feval(func, x_dif, t, y);  %vector-valued function in the changed x vector
    J(:, i) = (F_dif - F)/dx;          %approximation of Jacobian(finite differences)
    
end


end
