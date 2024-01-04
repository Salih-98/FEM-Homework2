%% Linear FEM - Homework 2 - Task 2
% --- Author: Salih Sahbegovic
% --- Date: 04.01.2024.
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

a = 4 + (X+Z)/10;
b = 0.5 + (W+Z)/20;
c = 1+ (X+Y)/20;
t = 1;

%% Define material properties

E = (45 + W)*10^3;
v = 0.15 + 0.01*Y;
lambda = (1.33 + 0.06*Z);
alfa0 = (1 + X/5)*10^(-6);

%%  Define loading
pmax = a * 0.997 * 9.81;
theta1 = 35;
theta2 = (40+Y);

%% Other
thetaRef = 25-W;

%% Define coordinates of the nodes
x1 = b+c ; y1 = 0;
x2 = b; y2 = a;
x3 = 0; y3 = 0;

%% Define area of the triangle
A = a*(b+c)/2;

%% Define B matrix
Bu = (1/(2*A))*[y2-y3 y3-y1 0;0 0 x1-x3; x3-x2 x1-x3 y3-y1];

%% Define stiffness tensor - plain strain state
C = (E / ((1+v)*(1-2*v))) * [1-v v 0; v 1-v 0; 0 0 (1-2*v)/2];

%% Define Kuu matrix
Kuu = A * t * Bu' * C * Bu;


%% Define Btheta 
Btheta = [y2-y3 y3-y1 0;x3-x2 x1-x3 y3-y1];

%% Define lambda matrix
lambda = lambda * eye(2,2);

%% Define global reduced conductivity matrix
Ktheta = Btheta' * lambda * Btheta * A *t;

%% Define reduced global coupling matrix
Kutheta = (A*t)/3 * Bu' * C* alfa0 * eye(3,3);