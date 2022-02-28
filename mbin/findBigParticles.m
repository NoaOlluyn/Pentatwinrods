function [X] = findBigParticles(X,sz_min)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

X = X(find(X(:,4)>sz_min),:);

X = X(find(X(:,5)>sz_min),:);

X = X(find(X(:,6)>sz_min),:);

end

