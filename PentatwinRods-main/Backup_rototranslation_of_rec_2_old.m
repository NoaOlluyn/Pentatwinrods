%%%% This is the script that works for rec_2


clc; clearvars; close all;

addpath(genpath([pwd '/mbin/']));
addpath(genpath([pwd '/Tomohawk/']));

resDir = [pwd '/results/'];

%% load image
 
I = read_rec('needle_rec_2.rec');

%% threshold image

thr = 0.25;
It  = 0*I;
It(I > thr) = 1;
figure(3000)
volshow(It,'Renderer','Isosurface');
%% get connected components

CC = bwconncomp(It,18);


%% find right connected components

truePixelID = [];
min_voxels  = 1000;
S = regionprops(CC,'Centroid');

j = 1;
for i=1:length(CC.PixelIdxList)
    voxelID = CC.PixelIdxList{i};
    if length(voxelID) > min_voxels
        truePixelID{j} = voxelID;
        Scenter(j,:) = S(i).Centroid;
        j = j+1;
    end
end

%% generate Images

for i=1:length(truePixelID)
    Itp{i} = 0*It;
    voxelID = truePixelID{i};
    Itp{i}(voxelID) = It(voxelID);
    
    % Itp{i}(1:2,1:2,:) = 1;
    % Itp{i}(1:2,:,1:2) = 1;
    % Itp{i}(:,1:2,1:2) = 1;
    
    figure(100+i);volshow(Itp{i},'Renderer','Isosurface');
    pause(0.001);
end

%% generate bi-Cone

n       = size(It);
zmin    = -0.5;
zmax    = 0.5;
z       = linspace(zmin,zmax,n(1));
x       = linspace(zmin,zmax,n(2));
y       = linspace(zmin,zmax,n(3));
[zz,xx,yy] = ndgrid(z,x,y);
xyz     = [xx(:)';zz(:)';yy(:)'];

height = 0.4;                   % height
mid_rad= 0.1;                   % radius of the middle circle
angle  = tan(mid_rad/height);   % cone angle

% generate bi-Cone
X = biCone(xyz,height,angle);
X = reshape(X,n);

figure(2);volshow(X);

%% find the transformation matrix for the transformed image

fprintf('**** registering the volume \n');

% get 'optimizer' and 'metric'
[optimizer, metric] = imregconfig('multimodal');
optimizer.InitialRadius     = 6e-3;
optimizer.Epsilon           = 1.5e-6;
optimizer.GrowthFactor      = 1.05;
optimizer.MaximumIterations = 1000;


for i=1:length(Itp)
    
    It1 = imgaussfilt3(imresize3(Itp{i},0.2),2);
    X1  = imgaussfilt3(imresize3(X,0.2),2);
    
    % get the transformation matrix
    tform_reg{i} = imregtform(It1,X1,'rigid',optimizer, metric,...
        'DisplayOptimization',0,'PyramidLevels',3);
    
end

fprintf('volume registration complete \n');

%%

for i=1:length(tform_reg)
    Ti = tform_reg{i}.T(1:3,1:3);
    Q(i,1:3) =  Scenter(i,:) - n/2;
    Q(i,4:6) = rotmat2eulang(inv(Ti),'XYZ')*(180/pi);
end

%% aligned shape to the center:

for i=1:length(tform_reg)
    Ti = tform_reg{i}.T;
    Si = Ti;
    Si(1:3,1:3) = inv(Ti(1:3,1:3));
    Si(4,1:3)   = Scenter(i,:) - n/2; % 
    tform_reg1 = affine3d(Si);
    sameAsInput = affineOutputView(size(X),tform_reg1,'BoundsStyle','CenterOutput');
    [I_tran{i},rbp] = imwarp(X,tform_reg1,'OutputView',sameAsInput);
end

%%
Ir = 0*X;
for i=1:length(I_tran)
    Ir = Ir + I_tran{i};
end
figure(2000);volshow(Ir,'Renderer','VolumeRendering');

%% save results

save([resDir 'Q.mat'],'Q');

str = [];
str.It = logical(It);
str.Ir = logical(Ir);
str.Q  = Q;
str.X  = logical(X);
str.tf = tform_reg;
str.BiCone.height = height;
str.BiCone.radius = mid_rad;
save([resDir 'full.mat'],'str');
