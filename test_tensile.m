function [l,b,e,index_l,index_b] = test_tensile(resolution_l,resolution_b,solid_x,solid_y,length,radius)
% for loop to find the best coefficient
interation = 100;
error = zeros(interation,interation);
n = 0;
value = 1000000;
index = 0;
index_l = 0;
index_b = 0;
for i = 1:interation
    for j = 1:interation
        c_l = i*resolution_l/interation;
        c_b = j*resolution_b+1;
        n = n + 1;
        error(i,j) = find_tensile_coef(c_l,c_b,solid_x,solid_y,length,radius);
        if (error(i,j)<value)
            value = error(i,j);
            index = n;
            index_l = c_l;
            index_b = c_b;
        end
    end
end

l = resolution_l/100:resolution_l/100:resolution_l;
b = 1+resolution_b:resolution_b:1+100*resolution_b;
e = error;

end

function residual = find_tensile_coef(c_l,c_b,solid_x,solid_y,length,radius)
C = [1 2];
node = [0 0 0;
        length 0 0];
    
E = 150000; % E: Young's modulus of steel [N/mˆ2]
r = radius;
A = pi*r^2 ; % beam area [mˆ2]
s_node = node(C(1),:);  % s_node : Start node
e_node = node(C(2),:);  % n_node : end node
v = e_node - s_node;      % v : Vector of the beam element
l = norm(v);              % l : Length of the beam

l_j = c_l*l;
l_s = (1-2*c_l)*l;
K1 = E*A/l_j*[1,-1;-1,1];
K1 = c_b*K1;
K2 = E*A/l_s*[1,-1;-1,1];

Ke = zeros(4,4);
Ke(1:2,1:2) = Ke(1:2,1:2) + K1;
Ke(2:3,2:3) = Ke(2:3,2:3) + K2;
Ke(3:4,3:4) = Ke(3:4,3:4) + K1;

u_load = [0 0 0 1]';
u = zeros(4,1);
f = zeros(4,1);
f = f - Ke*u_load;
u(2:3) = Ke(2:3,2:3)\f(2:3);
u = u + u_load;
f = Ke*u;


% interpolation functions
phi = zeros(2,1);
u_in = zeros(31,1);
x_co = zeros(31,1);
% first segment
for i = 0:10
x = i*l_j/10;
he = l_j;




phi(1) = 1-x/he;
phi(2) =  x/he;


u_in(i+1) = u(1:2)'*phi;
x_co(i+1) = x;
end
% second segment

for i = 1:10
x = i*l_s/10;
he = l_s;


phi(1) = 1-x/he;
phi(2) =  x/he;


u_in(i+11) = u(2:3)'*phi;
x_co(i+11) = x+l_j;

end

% third segment
for i = 1:10
x = i*l_j/10;
he = l_j;

phi(1) = 1-x/he;
phi(2) =  x/he;


u_in(i+21) = u(3:4)'*phi;
x_co(i+21) = x+l_s+l_j;
end
x_co(31) = x_co(31) - 10^-10;
% plot(x_co,u_in)
% hold on

u_solid = interp1q(solid_x,solid_y,x_co);
% plot(x_co,u_solid);
% hold off
residual = norm(u_solid-u_in);




end