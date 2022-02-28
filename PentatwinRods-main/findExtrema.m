function[v,A,B]=findExtrema(I)
%input a 3d image I, gives back the 2 extremal points for the projection
%along the X-axis, the the Y axis, and the Z axis

Ixy=squeeze(sum(I,3)); %Project the image along the 3 axis to generate 3, 2D images
Ixz=squeeze(sum(I,2));
Iyz=squeeze(sum(I,1));

extremaXY=regionprops(Ixy,'Extrema'); %Find 8 extremal points, this gaves a 
extremaXZ=regionprops(Ixz,'Extrema'); %vector with [top-left top-right right-top 
extremaYZ=regionprops(Iyz,'Extrema'); %right-bottom bottom-right bottom-left left-bottom left-top].

A=extremaXY(1); %assign starting values for the end points
B=extremaXY(2);

for i=1:length(extremaXY)-1
    for j=i+1:length(extremaXY)
        pA=extremaXY(i);
        pB=extremaXY(j);
        
        pC=extremaXZ(i);
        pD=extremaXZ(j);
        
        pE=extremaYZ(i);
        pF=extremaYZ(j);
        if norm(pA-pB)>norm(A-B)
            A=pA;
            B=pB;
        end
        
        if norm(pC-pD)>norm(A-B)
            A=pC;
            B=pD;
        end
        
        if norm(pE-pF)>norm(A-B)
            A=pE;
            B=pF;
        end
    end
end
v=A-B;