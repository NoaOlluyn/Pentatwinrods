function [transParam] = rotoTransFilter(I,X,tau)
% roto-translation filter
% I is an actual image (for example, raw reconstruction image)
% X is a shape (e.g., biCone)


%% set 'optimizer' and 'metric'

[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations          = 100;          
optimizer.GradientMagnitudeTolerance = 1e-4;
optimizer.MinimumStepLength          = 1e-4;

%% find stride lengths, etc

% get sizes
nI = size(I);
nX = size(X);

% stride for moving the filter
stride = 4;

% number of total filtering
nB = floor((nI-nX)/stride);

% size of shape
sz_shape = nnz(X);

%% main

transParam = [];

shapeId = 1;

% loop over all strides
for i=1:nB(1)
    for j=1:nB(2)
        for k=1:nB(3)
            
            % get local volume
            I_stride = I((i-1)*stride+(1:nX(1)),(j-1)*stride+(1:nX(2)),(k-1)*stride+(1:nX(3)));
            
            % set max-value to 1
            maxIstr = max(I_stride(:));
            if maxIstr > 0
                I_stride = I_stride/maxIstr;
            end
            
            % check if there is one whole shape in the sample
            if nnz(I_stride) > tau * sz_shape
                
                % find roto-translation
                tform_reg = imregtform(I_stride,X,'rigid',optimizer, metric,...
                    'DisplayOptimization',0,'PyramidLevels',3);
                
                T = tform_reg.T;
                T(4,1:3) = T(4,1:3) + [(i-1) (j-1) (k-1)]*stride;
                
                transParam.T{shapeId} = T;
                shapeId = shapeId + 1;
            end 
        end
    end
end


end