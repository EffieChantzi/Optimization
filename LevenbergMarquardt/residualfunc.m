%% Chantzi Efthymia - Optimization - Assignment 1 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function provides the vector-valued objective function of the    %
% given test model: y(t)=x1*e^(x2*t).                                   %
%                                                                       %
%                                                                       %
% %%%% Inputs %%%%                                                      %
% x: vector of unknown parameters                                       %
% t: vector of t data                                                   %
% y: vector of y data                                                   %
%                                                                       %
% %%%% Outputs %%%%                                                     %
% F: vector-valued function-> F(x)= x1*e^(x2*t) - y(t)                  %
% J: matrix of Jacobian. It is an optional output, in order to test     %
% the function of the algorithm in both requested cases, user-supplied  %
% and not user-supplied.                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [F, J] = residualfunc(x, t, y)

F = x(1)*exp(x(2)*t) - y;

if(nargout > 1)
   
    J = [exp(x(2)*t)  x(1)*t.*exp(x(2)*t)]; 

end

end
