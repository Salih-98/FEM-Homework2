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

%% Define coordinates
x1=a; y1=0;
x2=a+b; y2=0;
x3=0; y3=a+b;
x4=0; y4=a;
x5=a+b/2; y5=0;
x6=(a+b)*cosd(45); y6=(a+b)*sind(45);
x7=0; y7=a+b/2;
x8=a*cosd(45); y8=a*sind(45);
x = [x1 y1; x2 y2; x3 y3; x4 y4; x5 y5; x6 y6; x7 y7; x8 y8];

%% Shape functions - derivatives with respect to x
syms k1; syms k2;
dN11 = (1/4)*(1-k2)*(2*k1+k2); dN12 = (1/4)*(1-k1)*(2*k2+k1);
dN21 = (1/4)*(1-k2)*(2*k1-k2); dN22 = (1/4)*(1+k1)*(2*k2-k1);
dN31 = (1/4)*(1+k2)*(2*k1+k2); dN32 = (1/4)*(1+k1)*(2*k2+k1);
dN41 = (1/4)*(1+k2)*(2*k1-k2); dN42 = (1/4)*(1-k1)*(2*k2-k1);

dN51 = -k1*(1-k2); dN52 = (-1/2)*(1-k1^2);
dN61 = (1/2)*(1-k2^2); dN62 = -(1+k1)*k2;
dN71 = -k1*(1+k2); dN72 = (1/2)*(1-k1^2);
dN81 = (-1/2)*(1-k2^2); dN82 = -(1-k1)*k2;

dNk = [dN11 dN21 dN31 dN41 dN51 dN61 dN71 dN81; ...
    dN12 dN22 dN32 dN42 dN52 dN62 dN72 dN82];

%% Gradient of shape with functions with respect to ksi1, ksi2 for integration points
dNk1 = double(subs(dNk, [k1, k2], [-1/sqrt(3), -1/sqrt(3)]));
dNk2 = double(subs(dNk, [k1, k2], [1/sqrt(3),  -1/sqrt(3)]));

%% Jacobian matrix 
J1 = dNk1 * x;
J2 = dNk2 * x;

%% Inverse of jacobian matrix
invJ1 = inv(J1);
invJ2 = inv(J2);

%% Determinant of the jacobian matrix
detJ1 = det(J1);
detJ2 = det(J2);

%% Gradient of shape functions with respect to x
dNx1 = invJ1*dNk1;
dNx2 = invJ2*dNk2;

%% B Matrix 
B1 = [dNx1(1,1) 0 dNx1(1,2) 0 dNx1(1,3) 0 dNx1(1,4) 0 ...
      dNx1(1,5) 0 dNx1(1,6) 0 dNx1(1,7) 0 dNx1(1,8) 0 ;...
      0 dNx1(2,1) 0 dNx1(2,2) 0 dNx1(2,3) 0 dNx1(2,4)...
      0 dNx1(2,5) 0 dNx1(2,6) 0 dNx1(2,7) 0 dNx1(2,8);...
      dNx1(2,1) dNx1(1,1) dNx1(2,2) dNx1(1,2)... 
      dNx1(2,3) dNx1(1,3) dNx1(2,4) dNx1(1,4)...
       dNx1(2,5) dNx1(1,5) dNx1(2,6) dNx1(1,6)... 
      dNx1(2,7) dNx1(1,7) dNx1(2,8) dNx1(1,8)];

B2 = [dNx2(1,1) 0 dNx2(1,2) 0 dNx2(1,3) 0 dNx2(1,4) 0 ...
      dNx2(1,5) 0 dNx2(1,6) 0 dNx2(1,7) 0 dNx2(1,8) 0;...
      0 dNx2(2,1) 0 dNx2(2,2) 0 dNx2(2,3) 0 dNx2(2,4)...
      0 dNx2(2,5) 0 dNx2(2,6) 0 dNx2(2,7) 0 dNx2(2,8);...
      dNx2(2,1) dNx2(1,1) dNx2(2,2) dNx2(1,2)... 
      dNx2(2,3) dNx2(1,3) dNx2(2,4) dNx2(1,4)...
       dNx2(2,5) dNx2(1,5) dNx2(2,6) dNx2(1,6)... 
      dNx2(2,7) dNx2(1,7) dNx2(2,8) dNx2(1,8)];

%% Reduce the B matrix
B1 = getReducedMatrix(B1);
B2 = getReducedMatrix(B2);
B1 = applyConstraints(B1);
B2 = applyConstraints(B2);

%% Stiffness tensor C - plane stress
C = (E/(1-v^2))*[1 v 0; v 1 0; 0 0 (1-v)/2];

%% Stiffness matrix
K1 = 1*1*B1'*C*B1*t*detJ1;
K2 = 1*1*B2'*C*B2*t*detJ2;
K = 2*(K1+K2);

%% Define load vector
dX1 = dN22*x2 + dN32*x3;
dX2 = dN22*y2 + dN32*y3;
dGama = sqrt(dX1^2 + dX2^2);

N=[-1/2*(1-k2)*k2;(1-k2^2);1/2*(1+k2)*k2];

N1 = double(subs(N, k2, -sqrt(3/5)));
N2 = double(subs(N, k2, -0));
N3 = double(subs(N, k2, sqrt(3/5)));

dGama1 = double(subs(dGama, [k1, k2], [1, -sqrt(3/5)]));
dGama2 = double(subs(dGama, [k1, k2], [1, -0]));
dGama3 = double(subs(dGama, [k1, k2], [1, sqrt(3/5)]));

alfa11 = 5/9;
alfa22 = 8/9;
alfa33 = 5/9;

rLoc1 = t*(alfa11 * N1 * pc * dGama1);
rLoc2 = t*(alfa22 * N2 * pc * dGama2);
rLoc3 = t*(alfa33 * N3 * pc * dGama3);
rLoc = rLoc1 + rLoc2 + rLoc3;
rLoc = [0; rLoc(1);  0; rLoc(3); 0; rLoc(2)];
%% Transform local to global vector
alfa2 = 90;
alfa3 = 180;
alfa6 = 135;
T = [cosd(alfa2) sind(alfa2) 0 0 0 0;...
    -sind(alfa2) cosd(alfa2) 0 0 0 0; ...
    0 0 cosd(alfa6) sind(alfa6)  0 0; ...
    0 0 -sind(alfa6) cosd(alfa6) 0 0; ...
    0 0 0 0 cosd(alfa3) sind(alfa3); ...
    0 0 0 0 -sind(alfa3) cosd(alfa3)];
T (:, 5) = [];
T (:, 2) = [];
rGlobal = T' * rLoc;
rGlobal (1) = rGlobal(4) + rGlobal(1);
rGlobal (2) = rGlobal(2) + rGlobal(3);
rGlobal(4) = [];
rGlobal(3) = [];
Rglobal = [0; rGlobal(1); 0 ;rGlobal(2); 0];
%% Define displacements
u = inv(K) * Rglobal;

%% Define stresses
Sigma = C*B1*u;

function B1=getReducedMatrix(B1)
B1(:,13)=[];
B1(:,10)=[];
B1(:,7)=[];
B1(:,5)=[];
B1(:,4)=[];
B1(:,2)=[];
end

function B = applyConstraints(B)
B(:,1) = B(:,1) + B(:,4);
B(:,2) = B(:,2) + B(:,3);
B(:,5) = B(:,5) + B(:,8);
B(:,6) = B(:,6) + B(:,7);
B(:,9) = B(:,9) + B(:,10);
B(:,10)=[];B(:,8)=[];B(:,7)=[]; B(:,4)=[];B(:,3)=[];  
end


