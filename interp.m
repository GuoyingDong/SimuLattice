function modulus = interp(Xq,Yq,Zq,address)
length = [5 10 20 30];
diameter = [1 2 3];
angle = [0 15 30 45 60 75 90];
X = zeros(3,4,6);
Y = zeros(3,4,6);
Z = zeros(3,4,6);
V = zeros(3,4,6);
M = dlmread(address);
for i = 1:4
    for j = 1:3
        for k = 1:7
            X(j,i,k) = length(i);
            Y(j,i,k) = diameter(j);
            Z(j,i,k) = angle(k);
            for n = 1:84
                if M(n,1) == length(i)
                    if M(n,2) == diameter(j)
                        if M(n,3) == angle(k)
                            V(j,i,k) = M(n,4);
                        end
                    end
                end
            end
        end
    end
end
modulus = interp3(X,Y,Z,V,Xq,Yq,Zq,'spline');
end
            






            