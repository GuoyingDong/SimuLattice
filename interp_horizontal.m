function modulus = interp_horizontal(Xq,Yq,address)

length = [5 10 20 30];
diameter = [1 2 3];
X = zeros(3,4);
Y = zeros(3,4);
V = zeros(3,4);
M = dlmread(address);
for i = 1:4
    for j = 1:3
        X(j,i) = length(i);
        Y(j,i) = diameter(j);
        for n = 1:84
            if M(n,1) == length(i)
                if M(n,2) == diameter(j)
                    if M(n,3) == 0
                        V(j,i) = M(n,4);
                    end
                end
            end
        end
    end
end
modulus = interp2(X,Y,V,Xq,Yq,'spline');
end
