%% find the transformation matrix for 1 of the needles
% this is to test and to save time

fprintf('**** registering the volume \n');

% get 'optimizer' and 'metric'
[optimizer, metric] = imregconfig('multimodal');

optimizer.MaximumIterations = 10000;
% optimizer.InitialRadius     = 6e-3;
% optimizer.Epsilon           = 1.5e-6;
% optimizer.GrowthFactor      = 1.05;


for i=1:length(Itpp)
    transX       = (Scenter(i,:)*tau - n/2)';
    tform        = affine3d([eye(3) zeros(3,1);transX' 1]);
    X_translated = imwarp(X,tform); % get the translated bicone
    It1          = Itpp{i};         % imgaussfilt3(Itpp{i},2);
    X1           = X_translated;    % imgaussfilt3(X_translated,2);
    
    % get the transformation matrix
    tform_reg{i} = imregtform(It1,X1,'rigid',optimizer, metric,...
        'DisplayOptimization',0,'PyramidLevels',3);
    
end

fprintf('volume registration complete \n');
