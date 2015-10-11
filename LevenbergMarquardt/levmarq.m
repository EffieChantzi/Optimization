%% Chantzi Efthymia - Optimization - Assignment 1 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the function that implements the Levenberg-Marquardt algorithm     %
% for an arbitrary nonlinear least squares problem.                          %
%                                                                            %
%                                                                            %
% %%%% Inputs %%%%                                                           %
% func: defines the least-squares problem (vector-valued function)           %
% x0: vector of initial guess/solution                                       %
% t: vector of t data                                                        %
% y: vector of y data                                                        %                                                                                                                
% varargin: optional set of optimization parameters defined by user          %
% Otherwise, the following default values are set:                           %
%                                                                            % 
% %%%% Default Values %%%%                                                   %
% jacobianMode = 'off'->finite differences, no user-supplied gradient        %
% linear = 'on'->direction p by solving the linear least squares equation    %
% tol = 1e-7->tolerance parameter for error stopping criterion               %
% mu = 1e2->damping factor of Levenberg-Marquardt                            %
% maxIter = 350->maximum number of iterations performed                      %
% The user can change any of these by passing a cell structure 'options'     % 
% as the last input argument(varargin)                                       %
%                                                                            %
% %%%% Outputs %%%%                                                          %
% x: solution vector                                                         %
% resnorm: residual norm at the solution x                                   %
% residual: residual vector at the solution x                                %
% iterationHistory: optional cell structure output argument if the user      %
% wants to be informed of all the intermediate solutions and norm of         %
% gradient towards convergence                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x, resnorm, residual, iterationHistory] = levmarq(func, x0, t, y, varargin)


if(nargin < 4) %function aborted if the number of input is less than 
               %the four necessary: objective function, initial guess, t data, y data
   
    error('Invalid number of inputs');
        
else
    
    %default values of optimization parameters
    jacobianMode = 'off'; 
    linear = 'on';
    tol = 1e-7;
    mu = 1e2;
    maxIter = 350;
    displayIter = 'off';
   
    if(nargin == 5) %optimization option parameters are set/changed by the user
       
        options = varargin{:};
        if(length(options) > 12)
        
            error('Maximum number of optimization parameters exceeded.');

        else

            for i = 1 : 2 : (length(options) - 1) %check for valid/invalid options passed by user

                if(strcmpi(options{i}, 'jacobian') == 1)

                    if((strcmpi(options{i + 1}, 'on') == 1) || (strcmpi(options{i + 1}, 'off') == 1))
                        jacobianMode = options{i + 1};
                    else
                        error('Invalid Jacobian option');
                    end

                elseif(strcmpi(options{i}, 'linear') == 1)

                    if((strcmpi(options{i + 1}, 'off') == 1) || (strcmpi(options{i + 1}, 'on') == 1))
                        linear = options{i + 1};
                    else
                        error('Invalid least squares problem');
                    end

                elseif(strcmpi(options{i}, 'tolerance') == 1)

                    if(~isnumeric(options{i + 1}))
                        error('Invalid non-numeric tolerance option');
                    else
                        tol = options{i + 1};
                    end


                elseif(strcmpi(options{i}, 'dampingFactor') == 1)

                    if(~isnumeric(options{i + 1}))
                        error('Invalid non-numeric damping factor option');
                    else
                        if(options{i + 1} >= 0)
                            mu = options{i + 1};
                        else
                            error('Invalid negative damping factor')
                        end
                    end


                elseif(strcmpi(options{i}, 'iterations') == 1)

                    if(~isnumeric(options{i + 1}))
                        error('Invalid non-numeric maximum iterations option');
                    else
                        if(options{i + 1} >= 1)
                            maxIter = options{i + 1};
                        else
                            error('Invalid negative or zero iterations option');
                        end
                    end
                    
                elseif(strcmpi(options{i}, 'display') == 1)
                    
                    if(strcmpi(options{i + 1}, 'iter') == 1) || (strcmpi(options{i + 1}, 'off') == 1)
                        displayIter = options{i + 1};
                    else
                        error('Invalid iterative display option');
                    end

                else

                    error('Unidentified optimization parameter option');

                end

            end

        end     
        
    end
    
end


x = x0;
if(size(x, 2) > size(x, 1))  %the user may give the 
                             %initial guess either as a row or a column vector,                  
       
    x = x'; 
    
end

numOfParam = length(x);
I = eye(numOfParam, numOfParam);
iter = 0;
counter = 1;
xIterations(:, counter) = x;


if(strcmpi(jacobianMode, 'on') == 1)   %user supplies the gradient/Jacobian
    
    [r, J] = feval(func, x, t, y);    %user-supplied gradient/Jacobian 
        
    if(size(J, 1) == numOfParam)   %gradient transformed to Jacobian
        
        J = J';
            
    end
 
else             %numerical approximation of Jacobian
   
    r = feval(func, x, t, y);
    
    J = findJacobian(func, r, x, t, y);
    
end
f(counter) = 1/2*(r')*r; %f=1/2*r'*r
gradient_f = J'*r;       

if(strcmpi(linear, 'on') == 1)
    
    linearForm = 1;
    
else
    
    linearForm = 0;
    
end

condition(counter) = norm(gradient_f);
zeroMatrix = zeros(numOfParam, size(r, 2));
while((condition(counter) > tol) && (iter < maxIter))
    
    iter = iter + 1;
    counter = counter + 1;
    
    x_k = x;
    
    if(linearForm == 0)
        
        p_k = -((J'*J) + (mu*I))\gradient_f; %Levenberg-Marquardt formula
                                             %However, J'*J may be singular or close to singular.
        
    else
        
        A = [J ; (sqrt(mu).*I)]; 
        b = [-r ; zeroMatrix];
        p_k = A\b;               %direction p as a solution to the linear
                                 %least squares equation. Default option. 
        
    end
    x = x_k + p_k;

    
    if(strcmpi(jacobianMode, 'on') == 1)
 
        [r, J] = feval(func, x, t, y);
        f(counter) = 1/2*(r')*r;

    else
   
        r = feval(func, x, t, y);
        J = findJacobian(func, r, x, t, y);
    
    end
    
    
    f(counter) = 1/2*(r')*r;
    
    if(f(counter) < f(counter - 1))
       
        mu = mu/10;
       
    else
        
        x = x_k;
        
        if(strcmpi(jacobianMode, 'on') == 1)
 
            [r, J] = feval(func, x, t, y);

        else
   
            r = feval(func, x, t, y);
            J = findJacobian(func, r, x, t, y);
    
        end
        f(counter) = 1/2*(r')*r;
        mu = 10*mu;
        
    end
    gradient_f = J'*r;
    condition(counter) = norm(gradient_f);
    xIterations(:, counter) = x;
    
end

x = xIterations(:, end);
resnorm = 2*f(end);
residual = r;

if ((nargout > 3) && (strcmpi(displayIter, 'iter') == 1))
   
    iterationHistory{:, 1} = xIterations;
    iterationHistory{:, 2} = condition;
                     
end

end
