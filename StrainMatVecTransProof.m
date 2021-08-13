syms tx ty tz

% Generates a 6D strain tranformation matrix from 
% a 3D vector tranformation matrix
N = @(a)[a.*a ...
    a(:,[3 1 2]) .* a(:,[2 3 1]); ...
    a([3 1 2],:) .* a([2 3 1],:) .* 2 ...
    a([3 1 2],[3 1 2]) .* a([2 3 1],[2 3 1]) + a([3 1 2],[2 3 1]) .* a([2 3 1],[3 1 2])];

% Tranforms a strain matrix to a 6-vector representation
T = @(S)[S(1,1);S(2,2);S(3,3);S(2,3)+S(3,2);S(1,3)+S(3,1);S(1,2)+S(2,1)];

Ax = [1 0 0; 0 cos(tx) sin(tx); 0 -sin(tx) cos(tx)];
Ay = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
Az = [cos(tz) sin(tz) 0; -sin(tz) cos(tz) 0; 0 0 1];

S = sym('S%d%d', [3 3]);

A = sym('A%d%d', [3 3]);

simplify(T(A*S*A.') - N(A)*T(S))


