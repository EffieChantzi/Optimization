%% Chantzi Efthymia - Optimization - Assignment 3 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STAFF PLANNING                                                             % 
% This is the .m file that implements all the questions requested in         % 
% assignment 3. If run, the results are saved in a .txt file named           %
% output.txt.                                                                %                                                             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear all;
close all;

diary on
diary('output.txt'); %save output displayed results in output.txt


%% Task 1 
c = [2 2 2 3 4 3]';

A = [-1 0 0 0 0 -1
     -1 -1 0 0 0 0
     0 -1 -1 0 0 0
     0 0 -1 -1 0 0
     0 0 0 -1 -1 0
     0 0 0 0 -1 -1];

b = [-700 -500 -600 -300 -100 -50]';
lb = zeros(6, 1);
%x0 = zeros(6, 1);

fprintf('---- Task 1 ----');
fprintf('\n');
%options = optimoptions(@linprog, 'Algorithm', 'dual-simplex');
options = optimset('Algorithm', 'active-set', 'Largescale', 'off');
[x, fval, exitflag, output, lambda] = linprog(c, A, b, [], [], lb, [], [], options)
fprintf('\n');
fprintf('The minimum number of employees in each one of the six shifts is:\n');
y = [x(1) + x(6) ; x(1) + x(2) ; x(2) + x(3) ; x(3) + x(4) ; x(4) + x(5) ; x(5) + x(6)]

%% Task 2
b2 = [-700 -250 -600 -300 -100 -50]';
fprintf('---- Task 2 ----\n');
fprintf('---- # of staff shift2 has been reduced to half (250)..\n');
fprintf('\n');

%options = optimoptions(@linprog, 'Algorithm', 'dual-simplex');
[x, fval, exitflag, output, lambda] = linprog(c, A, b2, [], [], lb, [], [], options)
fprintf('\n');
fprintf('The minimum number of employees in each one of the six shifts is:\n');
y = [x(1) + x(6) ; x(1) + x(2) ; x(2) + x(3) ; x(3) + x(4) ; x(4) + x(5) ; x(5) + x(6)]

%% Task 3
fprintf('---- Task 3 ----\n');
fprintf('\n');

options = optimset('Largescale', 'on');
[x, fval, exitflag, output, lambda] = linprog(c, A, b, [], [], lb)%, [], [], options)
fprintf('The minimum number of employees in each one of the six shifts is:\n');
y = [x(1) + x(6) ; x(1) + x(2) ; x(2) + x(3) ; x(3) + x(4) ; x(4) + x(5) ; x(5) + x(6)]



diary off
