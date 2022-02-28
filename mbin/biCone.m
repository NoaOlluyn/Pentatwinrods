function [y] = biCone(x,s,alpha)

% axis coordinates for grid
xy= x(1:2,:);   % x- and y- coordinate
z = x(3,:);     % z-coordinate

% pre-processing
sxy = sqrt(sum(xy.^2,1));       % summing x- and y- coordinates
z1  = tan(alpha)*max(s-z,0);    % radius for upper-cone
z2  = tan(alpha)*max(s+z,0);    % radius for lower-cone

% upper-cone
y1 = zeros(1,size(x,2));
y1(sxy < z1) = 1;
y1(z<0) = 0;

% lower-cone
y2 = zeros(1,size(x,2));
y2(sxy < z2) = 1;
y2(z>0) = 0;

% join two cones
y = y1 + y2;
y = y(:);

end



