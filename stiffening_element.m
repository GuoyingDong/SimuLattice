function [u, force] = stiffening_element(strut_list, node_list, radius, ...
    youngs, poisson, bending_s, bending_l, tensile_s, tensile_l, load, support)
top = load;
fixdofs = support;
% *****Beam Element Solver*****
% Recommend units:
% mm, MPa, N
% =========input of the data===================
C = strut_list;
node = node_list;
% Now only consider homogeneous lattice structure
r = radius;
% =========material property===================
E = youngs; % E: Young's modulus of steel [N/mˆ2]
nu = poisson; % nu: Poisson's ratio
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
c_l = bending_l;  % c_l is the coefficient of the joint length
c_al = tensile_l;
c_b = bending_s;
c_t = tensile_s;
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
    
    %%%%%% determine the orientation of the strut %%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
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
   K_j = zeros(4,4,2);            % first 4*4 for Z, second 4*4 for Y
    K_j(:,:,1) = [K11(2,2)  K11(2,6)   -K11(2,2)  K11(2,6);
                  K11(2,6)  K11(6,6)   -K11(2,6)  K21(6,6);
                 -K11(2,2) -K11(2,6)    K11(2,2) -K11(2,6);
                  K11(2,6)  K21(6,6)   -K11(2,6)  K11(6,6)];

    K_j(:,:,2) = [K11(3,3)  K11(3,5)   -K11(3,3)  K11(3,5);
                  K11(3,5)  K11(5,5)   -K11(3,5)  K21(5,5);
                 -K11(3,3) -K11(3,5)    K11(3,3) -K11(3,5);
                  K11(3,5)  K21(5,5)   -K11(3,5)  K11(5,5)];
    K_j = c_b * K_j;
    % K_s is the stiffness matrix for the strut rotating around Z axis
    K_s = zeros(4,4,2);            % first 4*4 for Z, second 4*4 for Y
    K_s(:,:,1) = [K11_s(2,2)  K11_s(2,6)   -K11_s(2,2)  K11_s(2,6);
                  K11_s(2,6)  K11_s(6,6)   -K11_s(2,6)  K21_s(6,6);
                 -K11_s(2,2) -K11_s(2,6)    K11_s(2,2) -K11_s(2,6);
                  K11_s(2,6)  K21_s(6,6)   -K11_s(2,6)  K11_s(6,6)];

    K_s(:,:,2) = [K11_s(3,3)  K11_s(3,5)   -K11_s(3,3)  K11_s(3,5);
                  K11_s(3,5)  K11_s(5,5)   -K11_s(3,5)  K21_s(5,5);
                 -K11_s(3,3) -K11_s(3,5)    K11_s(3,3) -K11_s(3,5);
                  K11_s(3,5)  K21_s(5,5)   -K11_s(3,5)  K11_s(5,5)];
    
    K_c = zeros(4,4,2);
    for axis = 1:2  % 1 for Z axis, 2 for Y axis
        K1 = [K_j(1,1,axis) K_j(1,2,axis) 0             0;
              K_j(2,1,axis) K_j(2,2,axis) 0             0;
              0             0             K_j(3,3,axis) K_j(3,4,axis);
              0             0             K_j(4,3,axis) K_j(4,4,axis)];
        K2 = [K_j(1,3,axis) K_j(1,4,axis) 0             0;
              K_j(2,3,axis) K_j(2,4,axis) 0             0;
              0             0             K_j(3,1,axis) K_j(3,2,axis);
              0             0             K_j(4,1,axis) K_j(4,2,axis)];
        K3 = [K_j(3,1,axis) K_j(3,2,axis) 0             0;
              K_j(4,1,axis) K_j(4,2,axis) 0             0;
              0             0             K_j(1,3,axis) K_j(1,4,axis);
              0             0             K_j(2,3,axis) K_j(2,4,axis)];
        K4 = [K_j(3,3,axis) K_j(3,4,axis) 0             0;
              K_j(4,3,axis) K_j(4,4,axis) 0             0;
              0             0             K_j(1,1,axis) K_j(1,2,axis);
              0             0             K_j(2,1,axis) K_j(2,2,axis)];
        K4 = K4 + K_s(:,:,axis);
        K_c(:,:,axis) = K1+K2*(K4\(-K3));
    end
    
    Ke11 = zeros(6,6);
    Ke = zeros(12,12);
    
    set_x = [1 4 7 10];
    set_y = [3 5 9 11];
    set_z = [2 6 8 12];
    
    % increase stiffness for tension and torsion
    l_aj = l*c_al;
    l_as = l*(1-2*c_al);
    
    k_t = 1/(l_aj/(E*A*c_t) + l_as/(E*A) + l_aj/(E*A*c_t));
    G_t = 1/(l_j/(G*It*c_b) + l_s/(G*It) + l_j/(G*It*c_b));
    
%     Ke_x = [E*A/l 0      -E*A/l  0;
%             0     G*It/l  0     -G*It/l;
%            -E*A/l 0       E*A/l  0;
%             0    -G*It/l  0      G*It/l];
        
    Ke_x = [k_t 0      -k_t  0;
            0     G_t  0     -G_t;
           -k_t 0       k_t  0;
            0    -G_t  0      G_t];
    
    Ke(set_x,set_x) = Ke_x;
    Ke(set_z,set_z) = K_c(:,:,1);
    Ke(set_y,set_y) = K_c(:,:,2);
    
    % Rotation matrix
    l_11 = (e_node(1)-s_node(1))/l;
    l_21 = (e_node(2)-s_node(2))/l;
    l_31 = (e_node(3)-s_node(3))/l;
    D = sqrt(l_11*l_11 + l_21*l_21);
    if (D == 0)
        R = [0,0,1;0,1,0;-1,0,0];
    else
    l_12 = -l_21/D;
    l_22 = l_11/D;
    l_32 = 0;
    l_13 = -l_11*l_31/D;
    l_23 = -l_21*l_31/D;
    l_33 = D;
    R = [l_11 l_21 l_31;
        l_12 l_22 l_32;
        l_13 l_23 l_33];
    T = zeros(12,12);
    end
    T(1:3,1:3) = R;
    T(4:6,4:6) = R;
    T(7:9,7:9) = R;
    T(10:12,10:12) = R;
    
    % Rotate local stiffness matrix
    Ke = T' * Ke * T;
    
    % Assemble the global stiffness matrix
    edof = zeros(12,1);  % edof: Local node number to global node number
    for j = 1:6
        edof(j) = C(i,1)*6-6+j;
        edof(j+6) = C(i,2)*6-6+j;
    end
    for ii = 1:12
        for jj = 1:12
            count = count + 1;
            iK(count) = edof(ii);
            jK(count) = edof(jj);
            sK(count) = Ke(ii,jj);
        end
    end
end
DOF = N_node * 6;   % DOF: Total degrees of freedom
K = sparse(iK,jK,sK,DOF,DOF); % K: The global stiffness matrix
f = zeros(DOF,1);             % f: The load vector
u = zeros(DOF,1);             % u: The nodal displacement 
% ****Force load****
% f(load(:,1)) = load(:,2);     % Apply the load to the load vector
% ****Displacement load****
u_load = zeros(size(node_list,1)*6,1); % u_load: displacement load
u_load(top*6-3) = -1;                % top nodes with -1 z displacement
f = f - K*u_load;             % K*u_remain = f - K*u_load
activedofs = 1:DOF;
activedofs(fixdofs) = [];
u(activedofs) = K(activedofs,activedofs)\f(activedofs);
u = u + u_load;
f = K*u;
force = sum(f(top*6-3));
end




