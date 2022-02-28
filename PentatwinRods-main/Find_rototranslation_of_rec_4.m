%%%% This is the script that works for rec_4
% This script determines the positions and orientations of all particles in
% a sparse self-assembly, built from particles with axial symmetry.

clc; clearvars; close all;

addpath(genpath([pwd '/mbin/']));
addpath(genpath([pwd '/Tomohawk/']));

resDir = [pwd '/results/'];

tau  = 0.2;         % resizing factor
plotFig        = 1; %Set to 1 for images of the whole assembly
plotAllNeedles = 0; %Set to 1 for images of all nanoparticles seperately

%% Import data
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

 
I = read_rec('needle_rec_4.rec');

%% threshold image

thr         = 0.30;         % threshold value
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
min_voxels  = 1000;             % minimum voxels for conn comp
S = regionprops(CC,'Centroid'); % get centroids of these conn comp

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
%% Manual cropping the connected components where two or more are still present
%Carefull! if you change the tau and or the threshold value these will no
%longer be correct!

%First I crop the components with 2 or more needles in them.

problem{1}=IccSub{27};
problem{1}(40:70,10:48,20:50) = 0;
%figure(51);volshow(problem{1},'Renderer','Isosurface');

problem{2}=IccSub{27};
problem{2}(40:70,48:70,20:50) = 0;
%figure(52);volshow(problem{2},'Renderer','Isosurface');

problem{3}=IccSub{33};
problem{3}(1:70,10:37,20:50) = 0;
%figure(53);volshow(problem{3},'Renderer','Isosurface');

problem{4}=IccSub{33};
problem{4}(1:70,37:60,20:50) = 0;
%figure(54);volshow(problem{4},'Renderer','Isosurface');

problem{5}=IccSub{50};
problem{5}(20:60,10:70,20:50) = 0;
%figure(55);volshow(problem{5},'Renderer','Isosurface');

problem{6}=IccSub{50};
problem{6}(1:20,10:70,20:50) = 0;
%figure(56);volshow(problem{6},'Renderer','Isosurface');

problem{6}=IccSub{50};
problem{6}(1:20,10:70,20:50) = 0;
%figure(56);volshow(problem{6},'Renderer','Isosurface');

problem{7}=IccSub{51};
problem{7}(1:46,1:70,20:50) = 0;
%figure(57);volshow(problem{7},'Renderer','Isosurface');

problem{8}=IccSub{51};
problem{8}(46:70,1:70,20:50) = 0;
%figure(58);volshow(problem{8},'Renderer','Isosurface');

problem{9}=IccSub{67};
problem{9}(1:60,20:70,20:62) = 0;
%figure(59);volshow(problem{9},'Renderer','Isosurface');

problem{10}=IccSub{67};
problem{10}(1:60,1:19,20:62) = 0;
%figure(60);volshow(problem{10},'Renderer','Isosurface');

problem{11}=IccSub{73};
problem{11}(1:60,1:35,20:70) = 0;
problem{11}(1:28,1:60,20:70) = 0;
%figure(61);volshow(problem{11},'Renderer','Isosurface');

problem{12}=IccSub{73};
problem{12}(1:60,1:38,20:70) = 0;
problem{12}(29:50,1:60,20:70) = 0;
%figure(62);volshow(problem{12},'Renderer','Isosurface');

problem{13}=IccSub{73};
problem{13}(1:60,38:70,20:70) = 0;
problem{13}(22:50,1:60,20:70) = 0;
%figure(63);volshow(problem{13},'Renderer','Isosurface');

problem{14}=IccSub{73};
problem{14}(1:60,38:70,20:70) = 0;
problem{14}(1:21,1:60,20:70) = 0;
%figure(64);volshow(problem{14},'Renderer','Isosurface');

problem{15}=IccSub{81};
problem{15}(1:60,32:70,20:70) = 0;
%figure(65);volshow(problem{15},'Renderer','Isosurface');

problem{16}=IccSub{81};
problem{16}(1:60,1:31,20:70) = 0;
%figure(66);volshow(problem{16},'Renderer','Isosurface');

%Then we reassign the disconnected needles into our array of connected
%components. 

%% Insert the cropped CC into the list
%Then, reassign the disconnected needles into our array of connected
%components
IccSub{27}=problem{1};
IccSub{82}=problem{2}; 
IccSub{33}=problem{3};
IccSub{83}=problem{4}; 
IccSub{50}=problem{5};
IccSub{84}=problem{6}; 
IccSub{51}=problem{7};
IccSub{85}=problem{8}; 
IccSub{67}=problem{9}; 
IccSub{86}=problem{10}; 
IccSub{73}=problem{11}; 
IccSub{87}=problem{12}; 
IccSub{88}=problem{13}; 
IccSub{89}=problem{14}; 
IccSub{81}=problem{15}; 
IccSub{90}=problem{16}; 

%% Get Centroids
for i=1:length(IccSub)
stats = regionprops(IccSub{i}, 'Centroid');
C(i,:)=stats.Centroid;
end
Scenter=C;

%% Plot all IccSub in one image
Iy=0*IccSub{1};
for i = 1:length(IccSub)
    Iy = Iy + IccSub{i};
end
figure(501);volshow(Iy,'Renderer','Isosurface');

%% Plot all IccSub seperately
for i = 1:length(IccSub)
    if plotAllNeedles
        figure(100+i);volshow(IccSub{i},'Renderer','Isosurface');pause(0.01)
    end
end

%% generate bi-Cone
height = 0.18;                   % height
mid_rad= 0.04;                   % radius of the middle circle

Bicone = GenerateBicone(height, mid_rad, Isub); %Generates a 3D image of a bicone with parameters height and mid-rad

if plotFig, figure(2); volshow(Bicone); end
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

S_new = floor(Scenter);

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

save('Results_Rec4_28092021.mat','str');
%% write coördinates into excell
xlswrite('Scenter_rec_4.xlsx',Scenter)
%% 
xlswrite('thetai_rec_4.xlsx',thetai)
xlswrite('phii_rec_4.xlsx',phii)