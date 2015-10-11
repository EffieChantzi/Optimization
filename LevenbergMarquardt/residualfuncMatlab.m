%% Chantzi Efthymia - Optimization - Assignment 1 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function provides the vector-valued objective function of the    %
% given test model: y(t)=x1*e^(x2*t). It is written only for the built  %                                
% in function 'lsqnonlin', which does not accept y and t data as inputs,%                                                                      
% So, they have been defined inside.                                    %
%                                                                       % 
% %%% Inputs %%%%                                                        %
% x: vector of unknown parameters                                       %
%                                                                       %
% %%%% Outputs %%%%                                                     %
% F: vector-valued function-> F(x)= x1*e^(x2*t) - y(t)                  %
% J: matrix of Jacobian. It is an optional output, in order to test     %
% the function of the algorithm in both requested cases, user-supplied  %
% and not user-supplied.                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F, J] = residualfuncMatlab(x)

 y = [6.8 3.0 1.5 0.75 0.48 0.25 0.2 0.15]';
 t = [0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0]';

F = x(1)*exp(x(2)*t) - y;

if(nargout > 1)
   
    J = [exp(x(2)*t) ; x(1)*t.*exp(x(2)*t)]; 

end

end
