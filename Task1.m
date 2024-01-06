%% Linear FEM - Homework 2 - Task 1
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
Btheta = [y2-y3 y3-y1 y1-y2;x3-x2 x1-x3 x2-x1];

%% Define lambda matrix
lambda = lambda * eye(2,2);

%% Define global reduced conductivity matrix
Ktheta = Btheta' * lambda * Btheta * A *t;

%% Define reduced global coupling matrix
Kutheta = (A*t)/3 * Bu' * C* alfa0 * eye(3,3);

%% Define the consistent load vector
alfa = 90-atand(a/b);
L = sqrt(a^2 + b^2);
eqLoad2 = (pmax*L)/6;
Ru = [0; eqLoad2*cosd(alfa); -eqLoad2*sind(alfa)];

%% Static condensation
K = [Kuu -Kutheta; zeros(3,3) Ktheta];
rhs11 = Kutheta*thetaRef*ones(3,1);
rhs12 = Ru - rhs11;
Kn = getSortedMatrix(K);
rhs = [rhs12; 0] - Kn(1:4,5:6)*[theta1; theta2];
Knr = getReducedMatrix(Kn);
u = inv(Knr) * rhs;


function sortedMatrix = getSortedMatrix(K)
sortedMatrix = K;
sortedMatrix = [sortedMatrix(:, 1:3), sortedMatrix(:,6), sortedMatrix(:,4:5)];
sortedMatrix = [sortedMatrix(1:3,:); sortedMatrix(6,:); sortedMatrix(4:5,:)];
end

function reducedMatrix = getReducedMatrix(K)
reducedMatrix = K;
reducedMatrix(:, 5:6) = [];
reducedMatrix(5:6, :) = [];
end