% This script is to find the orientation of the particles that were too
% deformed by cupping artefacts and tresholding or rec_3, and should be run
% after Find_rototranslation_of_rec_3

%% Resize the full image
Isub = imresize3(I,tau);

%% Crop Isub around the deformed particle
Ic = Ir;
i=34;
x=S_new(i,2);
y=S_new(i,1);
z=S_new(i,3);
r=8;

Ic(1:x-r,:,:)=0;
Ic(:,1:y-r,:)=0;
Ic(:,:,1:z-r)=0;
Ic(x+r:90,:,:)=0;
Ic(:,y+r:90,:)=0;
Ic(:,:,z+r:90)=0;

if plotFig,figure(1);volshow(Ic,'Renderer','VolumeRendering');end

%%
i=54;
test=IccSub{i};
test(S_new(i,2),S_new(i,1),S_new(i,3))=1;
figure(100+i);volshow(test,'Renderer','Isosurface');pause(0.01)
