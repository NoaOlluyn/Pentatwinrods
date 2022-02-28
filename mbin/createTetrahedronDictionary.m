function [Dsub] = createTetrahedronDictionary(n,r,n_angles,tol)

tx = linspace(-pi,pi,n_angles);
ty = linspace(-pi,pi,n_angles);
tz = linspace(-pi,pi,n_angles);

t = 1;
D = zeros(n^3,n_angles^3,'logical');

fprintf('********* building dictionary ********* \n');
for i=1:n_angles
    for j=1:n_angles
        for k=1:n_angles
            X = tetrahedron(n,r,[tx(i) ty(j) tz(k)]);
            D(:,t) = X(:);
            t = t+1;
        end
    end
    fprintf('%d elements done | max=%d\n',i*(n_angles^2),n_angles^3);
end
fprintf('********* dictionary complete ********* \n');

if nargin < 4, tol = 1e-3;end

[Dsub,idx] = licols(single(D),tol);
Dsub       = logical(Dsub);
[nr,nc]    = size(Dsub);
%Dsub       = reshape(Dsub,n,n,n,nc);

%fprintf('Dictionary of size: %d x %d x %d x %d \n',size(Dsub));

end