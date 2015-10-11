%% Chantzi Efthymia - Optimization - Assignment 1 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main test file for function 'levmarq', as an alternative to call      %
% it directly in the command window. It requests the name of the        %
% objective function, provided that its .m file is placed in the        %
% folder. Results and deviations from the built in matlab function      %                                                                      
% 'lsqnonlin' are displayed in the command window and saved in a .txt   %                                                    
% file('output.txt').                                                   %                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all; 
close all;

diary on
diary('output.txt'); %save output displayed results in output.txt

format long g %format changed for displaying numbers with greater accuracy
fprintf('%s\n', date);
fprintf('---- Levenberg-Marquardt Algorithm ----\n');
str = input('Enter the name of your objective function: \n','s');
func = eval(['@' str]); %function handle
 
%y data (must be a column vector)
y = [6.8 3.0 1.5 0.75 0.48 0.25 0.2 0.15]'; 

%t data (must be a column vector)
t = [0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0]'; 
 
fprintf('\n\n');
fprintf('-------------------------------------------------------------------------------------------------\n');  


x = [200 -10]; %initial guess (either a column vector or row vector)
initial_guess = x';
fprintf('\n running...\n');
fprintf('--- levmarq Results, Initial guess: [%.2f %.2f] ----\n', initial_guess(1), initial_guess(2));
options = {'jacobian', 'off', 'display', 'iter'};
[x, resnorm, residual, hist] = levmarq(func, x, t, y, options)
history_x = hist{:, 1}'
history_norm = hist{:, 2}'
 
fprintf('-------------------------------------------------------------------------------------------------');
fprintf('\n\n\n');
 
fprintf('--- lsqnonlin Results, Initial guess: [%.2f %.2f] ----\n', initial_guess(1), initial_guess(2))
options = optimset('Display', 'iter-detailed', 'Algorithm', 'levenberg-marquardt');
[solMatlab, resnormMatlab, residualMatlab] = lsqnonlin(@residualfuncMatlab, x, [], [], options)

fprintf('-------------------------------------------------------------------------------------------------\n');
fprintf('\n\n');
fprintf('Comparison: ---levmarq--- & ---lsqnonlin---, Initial guess: [%.2f %.2f]\n', initial_guess(1), initial_guess(2));
fprintf('\n');

difference_x = abs(x - solMatlab)
difference_resnorm = abs(resnorm - resnormMatlab)
difference_residual = abs(residual - residualMatlab)
fprintf('%s\n', date);

diary off
