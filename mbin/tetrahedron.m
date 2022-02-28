function [X,meshFV] = tetrahedron(n,r,thetaV)


V1 = r*[-1;1;0;0];
V2 = r*[0;0;1;-1];
V3 = r*[-1/sqrt(2);-1/sqrt(2);1/sqrt(2);1/sqrt(2)];
F  = [1 2 3; 1 2 4; 2 3 4; 1 3 4];

V = [V1 V2 V3];

%%

if nargin < 3
    thetaV = [0 0 0];
end

thetaz = thetaV(3);
rotz   = [cos(thetaz)   -sin(thetaz) 0; ...
          sin(thetaz)    cos(thetaz) 0; ...
               0                  0  1];
thetay = thetaV(2);
roty   = [cos(thetay)    0  sin(thetay); ...
               0         1            0;...
          -sin(thetay)    0  cos(thetay)];

thetax = thetaV(1);
rotx   = [1           0                0;...
          0     cos(thetax) -sin(thetax); ...
          0     sin(thetax) cos(thetax)];

trans = [0, 0, 0];

rot   = rotz*roty*rotx;
tform = rigid3d(rot,trans);

ptCloudIn  = pointCloud(V);
ptCloudOut = pctransform(ptCloudIn,tform);

Vout = ptCloudOut.Location;
%%


meshFV.faces    = F;
meshFV.vertices = Vout;

x = linspace(-1,1,n);
y = linspace(-1,1,n);
z = linspace(-1,1,n);

X = VOXELISE(x,y,z,meshFV);

end