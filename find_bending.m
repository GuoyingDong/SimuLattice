function [x_co,u_in,u_solid] = find_bending(c_l,c_b,solid_x,solid_y,length,radius)
% *****Beam Element Solver*****
% Recommend units:
% mm, MPa, N
% =========input of the data===================
% C = [1 3;     % C: the connection of the node
%      1 2;
%      2 3;
%      2 4;
%      3 4];
% node = [0   0   0;   %node: node coordinate
%         10  0   0;
%         10  10  0;
%         20  0   0];
% fixdofs = [1:6 19:24];   % Boundary condition
% load = [14,-8;           % load: node id and force at that node
%         13,8];
C = [1 2];
node = [0 0 0;
        length 0 0];
top = 1;
bottom = 2;
% Now only consider homogeneous lattice structure
r = radius;
% =========material property===================
E = 150000; % E: Young's modulus of steel [N/mˆ2]
nu = 3/10; % nu: Poisson's ratio
G = E/2/(1+nu); % G: Shear modulus [N/mˆ2]
% =========calculate the beam property=========

A = pi*r^2 ; % beam area [mˆ2]
ASy = 9/10*A ; ASz = 9/10*A ; % effective area of shear
Iy = pi*r^4/4 ; Iz = Iy; % second moments of area [mˆ4]
It = 2*Iy; % torsion constant [mˆ4]
% Py = 0;
% Pz = 0;

N_element = size(C,1); % N: Total number of the elements
N_node    = size(node,1);
% Indexing vector
iK = zeros(144*N_element,1);   % iK: save i_th index
jK = zeros(144*N_element,1);   % jK: save j_th index
sK = zeros(144*N_element,1);   % sK: save the value in local stiffness matrix
count = 0;             % count: the index for iK jK sK

% Define the coefficient for the joint
% c_l = 0.15;  % c_l is the coefficient of the joint length
% c_b = 81.4;
c_m = 1;    % c_m is the coefficient of the moment of inertia

% Moment of inertia for the joint
Iz_j = Iz * c_m;
Iy_j = Iy * c_m;


% =========Get the stiffness matrix============
for i = 1:N_element
    s_node = node(C(i,1),:);  % s_node : Start node
    e_node = node(C(i,2),:);  % n_node : end node
    v = e_node - s_node;      % v : Vector of the beam element
    l = norm(v);              % l : Length of the beam
    
    
    % Length of the joint element
    l_j = c_l * l;
    % Length of the strut
    l_s = (1-2*c_l) * l;

    % Py and Pz for the joint
    Py_j = 12*E*Iz_j/(G*ASy*l_j^2) ; Pz_j = 12*E*Iy_j/(G*ASz*l_j^2) ; % Phi
    % Py and Pz for the strut
    Py = 12*E*Iz/(G*ASy*l_s^2) ; Pz = 12*E*Iy/(G*ASz*l_s^2) ; % Phi
    % Elements in the stiffness matrix for joint
    K11 = zeros(6,6);
    K21 = zeros(6,6);
    K11(2,2) = 12*E*Iz_j/(l_j^3*(1+Py_j)) ;
    K11(3,3) = 12*E*Iy_j/(l_j^3*(1+Pz_j)) ;
    K11(5,5) = (4+Pz_j)*E*Iy_j/(l_j*(1+Pz_j)) ;
    K11(6,6) = (4+Py_j)*E*Iz_j/(l_j*(1+Py_j)) ;
    K11(2,6) = 6*E*Iz_j/(l_j^2*(1+Py_j)) ;
    K11(3,5) = -6*E*Iy_j/(l_j^2*(1+Pz_j)) ;
    K21(6,6) = (2-Py_j)*E*Iz_j/(l_j*(1+Py_j));
    K21(5,5) = (2-Pz_j)*E*Iy_j/(l_j*(1+Pz_j)) ;
    
    % Elements in the stiffness matrix for strut
    K11_s = zeros(6,6);
    K21_s = zeros(6,6);
    K11_s(2,2) = 12*E*Iz/(l_s^3*(1+Py)) ;
    K11_s(3,3) = 12*E*Iy/(l_s^3*(1+Pz)) ;
    K11_s(5,5) = (4+Pz)*E*Iy/(l_s*(1+Pz)) ;
    K11_s(6,6) = (4+Py)*E*Iz/(l_s*(1+Py)) ;
    K11_s(2,6) = 6*E*Iz/(l_s^2*(1+Py)) ;
    K11_s(3,5) = -6*E*Iy/(l_s^2*(1+Pz)) ;
    K21_s(6,6) = (2-Py)*E*Iz/(l_s*(1+Py)) ;
    K21_s(5,5) = (2-Pz)*E*Iy/(l_s*(1+Pz)) ;

    % Calculate the stiffness matrix for bending
    % K_j is the stiffness matrix for the joint
    K_j = [K11(3,3)  K11(3,5)   -K11(3,3)  K11(3,5);
                  K11(3,5)  K11(5,5)   -K11(3,5)  K21(5,5);
                 -K11(3,3) -K11(3,5)    K11(3,3) -K11(3,5);
                  K11(3,5)  K21(5,5)   -K11(3,5)  K11(5,5)];
    K_j = c_b * K_j;

    % K_s is the stiffness matrix for the strut rotating around Z axis
    K_s = [K11_s(3,3)  K11_s(3,5)   -K11_s(3,3)  K11_s(3,5);
                  K11_s(3,5)  K11_s(5,5)   -K11_s(3,5)  K21_s(5,5);
                 -K11_s(3,3) -K11_s(3,5)    K11_s(3,3) -K11_s(3,5);
                  K11_s(3,5)  K21_s(5,5)   -K11_s(3,5)  K11_s(5,5)];
    
    Ke = zeros(8,8);
    Ke(1:4,1:4) = Ke(1:4,1:4) + K_j;
    Ke(3:6,3:6) = Ke(3:6,3:6) + K_s;
    Ke(5:8,5:8) = Ke(5:8,5:8) + K_j;
    
end
DOF = N_node * 4;   % DOF: Total degrees of freedom

f = zeros(DOF,1);             % f: The load vector
u = zeros(DOF,1);             % u: The nodal displacement 
% ****Force load****
% f(load(:,1)) = load(:,2);     % Apply the load to the load vector
% ****Displacement load****
u_load = zeros(DOF,1);
u_load(7) = 1;                % top nodes with -1 z displacement
f = f - Ke*u_load;             % K*u_remain = f - K*u_load
u(3:6) = Ke(3:6,3:6)\f(3:6);
u = u + u_load;
f = Ke*u;

% interpolation functions
phi = zeros(4,1);
u_in = zeros(31,1);
x_co = zeros(31,1);
% first segment
for i = 0:10
x = i*l_j/10;
he = l_j;
xa = 0;
x_bar = (x - xa)/he;
mu_e = 1 + Py_j;
gama_e = Py_j/12;

phi(1) = 1/mu_e * (mu_e - 12*gama_e*x_bar - (3 - 2*x_bar)*x_bar^2);
phi(2) = he/mu_e * ((1 - x_bar)^2 * x_bar + 6*gama_e*(1-x_bar)*x_bar);
phi(3) = 1/mu_e * ((3 - 2*x_bar)*x_bar^2 + 12*gama_e*x_bar);
phi(4) = he/mu_e * ((1 - x_bar) * x_bar^2 + 6*gama_e*(1-x_bar)*x_bar);

u_in(i+1) = u(1:4)'*phi;
x_co(i+1) = x;
end
% second segment

for i = 1:10
x = i*l_s/10+l_j;
he = l_s;
xa = l_j;
x_bar = (x - xa)/he;
mu_e = 1 + Py;
gama_e = Py/12;

phi(1) = 1/mu_e * (mu_e - 12*gama_e*x_bar - (3 - 2*x_bar)*x_bar^2);
phi(2) = he/mu_e * ((1 - x_bar)^2 * x_bar + 6*gama_e*(1-x_bar)*x_bar);
phi(3) = 1/mu_e * ((3 - 2*x_bar)*x_bar^2 + 12*gama_e*x_bar);
phi(4) = he/mu_e * ((1 - x_bar) * x_bar^2 + 6*gama_e*(1-x_bar)*x_bar);

u_in(i+11) = [0 0 u(5)+he*u(4)-u(3) 0]*phi+u(3)-(x-xa)*u(4);
x_co(i+11) = x;
end

% third segment
for i = 1:10
x = i*l_j/10+l_s+l_j;
he = l_j;
xa = l_s+l_j;
x_bar = (x - xa)/he;
mu_e = 1 + Py_j;
gama_e = Py_j/12;

phi(1) = 1/mu_e * (mu_e - 12*gama_e*x_bar - (3 - 2*x_bar)*x_bar^2);
phi(2) = he/mu_e * ((1 - x_bar)^2 * x_bar + 6*gama_e*(1-x_bar)*x_bar);
phi(3) = 1/mu_e * ((3 - 2*x_bar)*x_bar^2 + 12*gama_e*x_bar);
phi(4) = he/mu_e * ((1 - x_bar) * x_bar^2 + 6*gama_e*(1-x_bar)*x_bar);

u_in(i+21) = [0 0 u(7)+he*u(6)-u(5) -u(6)]*phi+u(5)-(x-xa)*u(6);
x_co(i+21) = x;
end
x_co(31) = x_co(31) - 10^-10;
% plot(x_co,u_in)
% hold on
u_solid = interp1q(solid_x,solid_y,x_co);
% plot(x_co,u_solid);
% hold off
residual = norm(u_solid-u_in);
end







