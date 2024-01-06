%% Linear FEM - Homework 2 - Task 2
% --- Author: Salih Sahbegovic
% --- Date: 05.01.2024.
% --- Subject: Finite Element Methods in Linear Structural Mechanics
% --- Semester: Winter Semester 2023/2024

%% Defining data given in the task formulation 

addpath('Functions\')
% In this section we first define W, X, Y, Z based on inputed
% immatriculation number, and then based on that value we define rest of
% given parameters a,b, c,d, A, E, pmax, prescribed values of F and u23, as
% well as spring stiffness 

% -- TODO: This is really rigid coding, just in order to perform
% calculations. Add more flexibility to the code later. 

% Extract W, X, Y, Z from the immatrikulation number
% TODO: Make it so the user can input it during execution, later.
iNum = '108022249956';
[W, X, Y, Z] = getWXYZ(iNum);

%% Define geometry 
a = 1 + (W+X)/10;
b = 5 + (Y+Z)/10;
t = 0.03 + W/300;

%% Define material properties
E = (18 + 0.2*X + 0.1*Y) * 10^4;
v = 0.23 + 0.01*Z;

%% Define loading 
pc = 5 + 0.3*W + 0.1*Y;