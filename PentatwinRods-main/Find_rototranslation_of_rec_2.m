%%%% This is the script that has the correct parameters for rec_2
% This script determines the positions and orientations of all particles in
% a sparse self-assembly, built from particles with axial symmetry.

clc; clearvars; close all;

addpath(genpath([pwd '/mbin/']));
addpath(genpath([pwd '/Tomohawk/']));

resDir = [pwd '/results/'];

tau  = 1;         % resizing factor
plotFig        = 1; %Set to 1 for images of the whole assembly
plotAllNeedles = 0; %Set to 1 for images of all nanoparticles seperately
%% Import data
%doesn't work
I=str.Image;
It=str.ThresholdedImage;

for i=1:length(str.ConnectedComponents)
IccSub{i} = str.ConnectedComponents{i};
A{i} = str.ExtremalPointsA{i};
B{i} = str.ExtremalPointsB{i};
v{i} = str.CentralAxisV{i};
thetai{i} = str.Theta{i};
phii{i} = str.phi{i};
rotmatMP{i} = str.RotMat{i};
end

Scenter = str.CentralPositions;
Bicone = str.X;
height = str.BiCone.height;
mid_rad = str.BiCone.radius;
Ir = str.Reconstruction;
n=str.ResizedSize;
%% load image through read_rec function of TomoHawk

 
I = read_rec('needle_rec_2.rec');

%% threshold image

thr         = 0.25;         % threshold value
It          = 0*I;          % initialize the thresholded image
It(I > thr) = 1;            % set value above threshold to 1.

% plot Figure
if plotFig,figure(1);volshow(It,'Renderer','Isosurface');end
%% get connected components


CONN = 18;                      % number of neighoring voxels to look for
CC   = bwconncomp(It,CONN);     % finds the connected components


%% find right connected components
% we only need ones that have size above min_voxels
% min_voxels may need to be adjusted manually to achieve the actual amount
% of NP's present

truePixelID = [];               % initialize the list (for connected comp)
min_voxels  = 3000;             % minimum voxels for conn comp
S = regionprops(CC,'Centroid'); % get central positions of these conn comp

% loop over all possible conn comp to find the right ones
j = 1;

for i=1:length(CC.PixelIdxList)
    
    voxelID = CC.PixelIdxList{i};
    
    if length(voxelID) > min_voxels     % condition 
        truePixelID{j} = voxelID;       % get voxel IDs
        Scenter(j,:)   = S(i).Centroid; % get their centroid
        j = j+1;
    end
end

%% Resize the images of the seperate particles, and generate images
%This greatly speeds up the calculation time


Isub = imresize3(It,tau);    % resize the image
n=size(Isub);

for i=1:length(truePixelID)
    
    Icc          = 0*It;            % initialize needle image
    voxelID      = truePixelID{i};  % get voxel IDs
    Icc(voxelID) = It(voxelID);     % set voxel IDs to true ones
    IccSub{i}      = imresize3(Icc,tau);  % resize needle
    
    if plotAllNeedles
        figure(100+i);volshow(IccSub{i},'Renderer','Isosurface');pause(0.01)
    end
    
    
end

% plot one needle (with axis)
if plotFig
    demonstration=IccSub{1};
    demonstration(1:3,1:3,:)=1;
    demonstration(1:3,:,1:3)=1;
    demonstration(:,1:3,1:3)=1;
    figure(101);volshow(demonstration,'Renderer','Isosurface');pause(0.01);
end
%% Plot all IccSub in one image
Iy=0*IccSub{1};
for i = 1:length(IccSub)
    Iy = Iy + IccSub{i};
end
figure(501);volshow(Iy,'Renderer','Isosurface');
%% generate bi-Cone
height = 0.4;                   % height
mid_rad= 0.1;                   % radius of the middle circle

Bicone = GenerateBicone(height, mid_rad, Isub); %Generates a 3D image of a bicone with parameters height and mid-rad

if plotFig, figure(2); volshow(Bicone); pause(0.001); end
%% Find the central axis of each rod
for i=1:length(IccSub)
    [A{i},B{i}]=findExtrema3d(IccSub{i});
    v{i}=A{i}-B{i};
end
%% calculate the angles

%loop over all particles
for i=1:length(v)
    %get the right coördinates and vectors
    vi=v{i};
    r=norm(vi);
    x=vi(1);
    y=vi(2);
    z=vi(3);
    
    %calculate the spherical angles
    theta=-atan(y/x);   %is the angle between the XY projection of v and the X-axis
    phi=-acos(z/r);     %is the angle between v and the Z-axis
    
    %store angles
    thetai{i}=theta;    
    phii{i}=phi;
    
%make SO(3) rotationmatrices
Sy=[cos(phi) 0 -sin(phi); 0 1 0; sin(phi) 0 cos(phi)];
Sz=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1];

%generate the final transformation matrix
rotmatMP{i}=Sy*Sz;
end

%% plot the roto-translated versions of biCones

S_new = floor(Scenter * tau);

% loop over all the needles
for i=1:length(rotmatMP)
    Ti              = rotmatMP{i};      % original transformation matrix
    Si              = Ti;               % new roto-translation matrix
    Si(4,1:3)       = S_new(i,:) - n/2; % set the translation
    Si(1:3,4)=0;                        % set the other elements
    Si(4,4)=1;                          
    tform_reg1      = affine3d(Si);     % use affine3d
    
    sameAsInput     = affineOutputView(size(Bicone),tform_reg1,...
                        'BoundsStyle','CenterOutput');  % output size
    
    % transform the bicones using the transofmation matrices
    [I_tran{i},rbp] = imwarp(Bicone,tform_reg1,'OutputView',sameAsInput);
end

% generate the total image
Ir = 0*Bicone;
for i = 1:length(I_tran)
    Ir = Ir + I_tran{i};
end

% plot total image
if plotFig
    figure(2000);volshow(Ir,'Renderer','VolumeRendering');
end
%% Plot all reconstructed needles individually
for i=1:length(I_tran)
    if plotAllNeedles
        figure(3000+i);volshow(I_tran{i},'Renderer','VolumeRendering');pause(0.01)
    end
end
%% save results
str = [];
str.Image = single(I);
str.ThresholdedImage = logical(It);
for i=1:length(I_tran)
str.ConnectedComponents{i} = logical(IccSub{i});
str.ExtremalPointsA{i} = double(A{i});
str.ExtremalPointsB{i} = double(B{i});
str.CentralAxisV{i} = double(v{i});
str.Theta{i} = double(thetai{i});
str.phi{i} = double(phii{i});
str.RotMat{i} = single(rotmatMP{i});
end
str.CentralPositions = double(Scenter);
str.X  = logical(Bicone);
str.BiCone.height = double(height);
str.BiCone.radius = double(mid_rad);
str.Reconstruction = logical(Ir);
str.ResizedSize = double(n);

save('Results_Rec1_23092021.mat','str');
%% write coördinates into excell
xlswrite('Scenter_rec_2.xlsx',Scenter)
xlswrite('thetai_rec_2.xlsx',thetai)
xlswrite('phii_rec_2.xlsx',phii)