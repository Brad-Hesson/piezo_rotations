syms d15 d22 d31 d33 tx ty tz real;

N = @(a)[a.*a ...
    a(:,[3 1 2]) .* a(:,[2 3 1]); ...
    a([3 1 2],:) .* a([2 3 1],:) .* 2 ...
    a([3 1 2],[3 1 2]) .* a([2 3 1],[2 3 1]) + a([3 1 2],[2 3 1]) .* a([2 3 1],[3 1 2])];

d = [0 0 0 0 d15 -2*d22;-d22 d22 0 d15 0 0;d31 d31 d33 0 0 0];

Ax = [1 0 0; 0 cos(tx) sin(tx); 0 -sin(tx) cos(tx)];
Ay = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
Az = [cos(tz) sin(tz) 0; -sin(tz) cos(tz) 0; 0 0 1];

dp = Ax*d*N(Ax).';

dparr = sqrt(dp(2,1)^2 + dp(2,3)^2);

fplot(dparr, [0 2*pi])

dp21 = expand(dp(2,1));
dp22 = expand(dp(2,2));
dp23 = expand(dp(2,3));
dp24 = expand(dp(2,4));
dp25 = expand(dp(2,5));
dp26 = expand(dp(2,6));


