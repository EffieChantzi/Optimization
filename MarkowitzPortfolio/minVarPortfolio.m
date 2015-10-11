%% Chantzi Efthymia - Optimization - Assignment 2 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PORTFOLIO OPTIMIZATION                                                     % 
% This is the .m file that implements all the questions requested in         % 
% assignment 2. If run, the results are saved in a .txt file named           %
% output.txt. The figure of efficient frontiers is also saved as a .png file %                                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

diary on
diary('output.txt'); %save output displayed results in output.txt

fprintf('%s\n', date);
fprintf('---- Portfolio Optimization ----\n');
fprintf('\n running...\n');

%%%%%%% general global parameters needed for all questions %%%%%%%

%number of assets
numOfAssets = 5;

%vector of yearly expected rates of returns of assets according to Table 1
annualExpectedRate = 10^(-2)*[13 5.3 10.5 5 12.6];

%symmetric matrix of covariances of assets according to Table 1
covarianceMatrix = 10^(-2)*[4.01 -1.19 0.6 0.74 -0.21;
                            -1.19 1.12 0.21 0.54 0.55;
                            0.6 0.21 3.04 0.77 0.29;
                            0.74 0.54 0.77 3.74 -1.04;
                            -0.21 0.55 0.29 -1.04 3.8];
                        
%change default options, so as to use a medium-scale active set algorithm
%as requested 
oldOptions = optimset('quadprog');
options = optimset(oldOptions, 'LargeScale', 'off');


%% Question 2 %%       

fprintf('---- Question 2, Results ----');
fprintf('\n\n');

%vector e for the definition of the constraint sum{i=1..n}w_i = 1
e = ones(numOfAssets, 1);                        

%parts of the reformulated system
LagragianRightSide = zeros(size(e));
constraintRightSide = 1;
b = [LagragianRightSide ; constraintRightSide];

%formulated system by the KKT conditiions
A = [covarianceMatrix -e ; e' 0];

%solving by backslash
w_lambda_sol = A\b

%variance 
variance_question2 = (1/2)*(w_lambda_sol((1 : 5), 1)'*covarianceMatrix*w_lambda_sol((1 : 5), 1)) %1/2w'Sw

%checks if the constraint is satisfied
constraintSatisfied = (w_lambda_sol(1 : 5, 1)'*e);

if((constraintSatisfied - 1) == 0)

    fprintf('The linear constraint is satisfied.\n');

end

%in order to validate the result, solve the same problem with quadprog
%not asked in the assignment, but it is done just for validation of the
%result
Aeq = [1 1 1 1 1];
beq = 1;
[w_quadprog, variance_quadprog] = quadprog(covarianceMatrix, [], [], [], Aeq, beq, [], [], [], options)



%% Question 3 %%
fprintf('\n\n');
fprintf('---- Question 3, Results ----');
fprintf('\n\n');


%% question 3b
p = 0.1;
fprintf('--------- p = %.2f ---------\n', p);

%short selling allowed                                   
Aeq = [annualExpectedRate ; e'];
beq = [p  1];

fprintf('---short selling---\n');
fprintf('weights, variance, expected rate\n');
[w, variance, exitflag] = quadprog(covarianceMatrix, [], [], [], Aeq, beq, [], [], [], options)
expectedReturn = annualExpectedRate*w


%short selling not allowed
lb = zeros(numOfAssets, 1); %lower bound, weights non-negative
fprintf('---no short selling---\n');
fprintf('weights, variance, expected rate\n');
[w, variance, exitflag] = quadprog(covarianceMatrix, [], [], [], Aeq, beq, lb, [], [], options)
expectedReturn = annualExpectedRate*w


%% question 3c
p = 0.2;
beq = [p 1];
fprintf('--------- p = %.2f ---------\n', p);

%short selling allowed                                   
fprintf('---short selling---\n');
fprintf('weights, variance, expected rate\n');
[w, variance, exitflag] = quadprog(covarianceMatrix, [], [], [], Aeq, beq, [], [], [], options)
expectedReturn = annualExpectedRate*w

%short selling not allowed
fprintf('---no short selling---\n');
fprintf('weights, variance, expected rate\n');
[w, variance, exitflag] = quadprog(covarianceMatrix, [], [], [], Aeq, beq, lb, [], [], options)
expectedReturn = annualExpectedRate*w


%% Question 4 %%
fprintf('\n\n');
fprintf('---- Question 4, Results ----');
fprintf('\n\n');

a = [0.05 : 0.05 : 1];

%equality constraint
Aeq = e';
beq = 1;

%weights for short selling allowed
weightsS = zeros(numOfAssets, length(a));

%matrix for storing variance and expected rate of return 
%of efficient frontier for plotting
efficientPairsS = zeros(2, length(a));  

fprintf('Short selling allowed.\n');
fprintf('\n');
for i = 1 : length(a)
     
    %inserting parameter a in the problem
    newCov = a(i)*covarianceMatrix;
    linearTerm = -(1 - a(i))*annualExpectedRate;
    
    [w, fval] = quadprog(newCov, linearTerm, [], [], Aeq, beq, [], [], [], options);
    
    %evaluation of variance because fval from above contains the value of
    %the reformulated problem described in question 4
    variancePort = (1/2)*w'*covarianceMatrix*w;
    
    %expected return of portfolio
    rPort = annualExpectedRate*w;
    
    %weights
    weightsS(:, i) = w;
    
    %efficient pairs
    efficientPairsS(1, i) = variancePort;
    efficientPairsS(2, i) = rPort;
     
end


%the same procedure as above applied for short selling not allowed

%weights for short selling not allowed
weights = zeros(numOfAssets, length(a));

%matrix for storing variance and expected rate of return 
%of efficient frontier for plotting
efficientPairs = zeros(2, length(a));  

fprintf('Short selling not allowed.\n');
fprintf('\n');
for i = 1 : length(a)

    %inserting parameter a in the problem
    newCov = a(i)*covarianceMatrix;
    linearTerm = -(1 - a(i))*annualExpectedRate;
    
    [w, fval] = quadprog(newCov, linearTerm, [], [], Aeq, beq, lb, [], [], options);
    
    %evaluation of variance because fval from above contains the value of
    %the reformulated problem described in question 4
    variancePort = (1/2)*w'*covarianceMatrix*w;
    
    %expected return of portfolio
    rPort = annualExpectedRate*w;
    
    %weights
    weights(:, i) = w;
    
    %efficient pairs
    efficientPairs(1, i) = variancePort;
    efficientPairs(2, i) = rPort;
    
end

%plot efficient pairs when short selling allowed and not allowed
figure();
plot(efficientPairsS(1, :), efficientPairsS(2, :), 'm*-', 'linewidth', 1.5);
hold on;
plot(efficientPairs(1, :), efficientPairs(2, :), '-.oc', 'linewidth', 1.5);
xlim([(min(efficientPairsS(1, :)) - 1) (max(efficientPairsS(1, :)) + 1)]);
ylim([(min(efficientPairsS(2, :)) - 1) (max(efficientPairsS(2, :)) + 1)])
xlabel('\sigma^2 (variance)', 'fontweight', 'bold');
ylabel('\mu (expected return)', 'fontweight', 'bold');
title('Efficient Frontier', 'fontweight', 'bold');
legend('short selling allowed', 'short selling not allowed', 'Location', 'best');
grid on;
hold off;
saveas(gcf, 'efficientPairsPortfolio.png');

fprintf('%s\n', date);
diary off
