syms tx ty tz real;

N = @(a)[a.*a ...
    a(:,[3 1 2]) .* a(:,[2 3 1]); ...
    a([3 1 2],:) .* a([2 3 1],:) .* 2 ...
    a([3 1 2],[3 1 2]) .* a([2 3 1],[2 3 1]) + a([3 1 2],[2 3 1]) .* a([2 3 1],[3 1 2])];

d = sym('d%d%d', [3 6]);

Ax = [1 0 0; 0 cos(tx) sin(tx); 0 -sin(tx) cos(tx)];
Ay = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
Az = [cos(tz) sin(tz) 0; -sin(tz) cos(tz) 0; 0 0 1];

Mx = [-1 0 0;0 1 0;0 0 1];
My = [1 0 0;0 -1 0;0 0 1];
Mz = [1 0 0;0 1 0;0 0 -1];

% the list of symmetries to be tested
dp = [];
%dp = [dp {subs(Ax*d*N(Ax).', tx, sym(2*pi/2))}];
%dp = [dp {subs(Ay*d*N(Ay).', ty, sym(2*pi/2))}];
dp = [dp {subs(Az*d*N(Az).', tz, sym(2*pi/3))}];
dp = [dp {Mx*d*N(Mx).'}];
%dp = [dp {My*d*N(My).'}];
%dp = [dp {Mz*d*N(Mz).'}];

% building the matrix to be reduced
mat = [];
for i = 1:length(dp)
    dp_t = dp(i);
    mat = [mat;sym_mat(d, dp_t{1})];
end

rmat = rref(mat);
eqs = nonzeros(rmat*d(:));
eqs_s = eqs;

% applying the equations to the d matrix
dout = d;
for i = 1:length(eqs)
    if length(symvar(eqs(i))) == 1
        dout = subs(dout, eqs(i), 0);
    else
        vars = symvar(eqs(i));
        var = vars(2);
        dout = subs(dout, var, solve(eqs(i), var));
        eqs = subs(eqs, var, solve(eqs(i), var));
    end
end
eqs = eqs_s;

% getting rid of fractional coefficients
for i = 1:length(dout)
    co = coeffs(dout(i));
    if abs(co) < 1
        var = symvar(dout(i));
        dout = subs(dout, var, var/co);
    end
end
% making sure the main diagonal is clean
for i = 1:3
    co = coeffs(dout(i,i));
    if co ~= 0
        dout = subs(dout, symvar(dout(i,i)), d(i,i)/co);
    end
end

dout

%--------------------------------------------------------------------------

% builds the reducing matrix from the before and after of a transformation
function mat = sym_mat(d_orig, d_trans)
    mat = [];
    for i = 1:18
        row = [];
        for j = 1:18
            % Gets the coefficient of "d_o(j)" in "d_t(i)-d_o(i)" and appends it to row.
            row = [row limit((d_trans(i)-d_orig(i))/d_orig(j), d_orig(j), Inf)];
        end
        mat = [mat; row];
    end
end