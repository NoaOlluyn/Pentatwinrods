function[A,B]=findExtrema3d(I)
%input a 3d image I, gives back the 2 extremal points for the projection
%along the X-axis, the the Y axis, and the Z axis

CH=regionprops3(I,'ConvexHull');  
voxelList=CH.ConvexHull{1};
A=squeeze(voxelList(1,:)); 
B=squeeze(voxelList(2,:));

for i=1:size(voxelList,1)-1
    for j=i+1:size(voxelList,1)
        pA=squeeze(voxelList(i,:));
        pB=squeeze(voxelList(j,:));
        if norm(pA-pB)>norm(A-B)
            A=pA;
            B=pB;
        end
    end
end
%A
%B
%distance=norm(A-B)