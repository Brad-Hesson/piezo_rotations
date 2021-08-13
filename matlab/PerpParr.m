syms d15 d22 d31 d33 tx ty tz x u v real;

N = @(a)[a.*a ...
    a(:,[3 1 2]) .* a(:,[2 3 1]); ...
    a([3 1 2],:) .* a([2 3 1],:) .* 2 ...
    a([3 1 2],[3 1 2]) .* a([2 3 1],[2 3 1]) + a([3 1 2],[2 3 1]) .* a([2 3 1],[3 1 2])];

d = [0 0 0 0 d15 -2*d22;-d22 d22 0 d15 0 0;d31 d31 d33 0 0 0];
d = subs(d, d15, sym('69.2e-12'));
d = subs(d, d22, sym('20.8e-12'));
d = subs(d, d31, sym('-0.85e-12'));
d = subs(d, d33, sym('6e-12'));

Eparr = @(d,r)dot(d(:,1:3).'*r,r./norm(r));
Eperp = @(d,r)norm(cross(d(:,1:3).'*r,r./norm(r)));
Sparr = @(d,r)dot(d(:,4:6).'*r,r./norm(r));
Sperp = @(d,r)norm(cross(d(:,4:6).'*r,r./norm(r)));

Ax = [1 0 0; 0 cos(tx) sin(tx); 0 -sin(tx) cos(tx)];
Ay = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
Az = [cos(tz) sin(tz) 0; -sin(tz) cos(tz) 0; 0 0 1];

r = Az.'*Ay.'*[1;0;0];
s = Ay*Az*d*N(Az).'*N(Ay).';

hold on
f = Eperp(subs(Ax*d*N(Ax).', tx, x),[0;1;0]);
fplot(f, [0 2*pi])
f = Sperp(subs(Ax*d*N(Ax).', tx, x),[0;1;0]);
fplot(f, [0 2*pi])
xline(128/180*pi)
hold off


figure()
f = Eparr(s,[1;0;0])*r;
f = subs(f, [ty tz], [u v]);
fsurf(f(1),f(2),f(3),[-pi/2 pi/2 -pi pi],'FaceAlpha',0.5);
%f = Eparr(s,[1;0;0]);
%fsurf(@(u,v)u,@(u,v)v,f,[0 2*pi/3 -pi/2 pi/2]);
title('Eparr') 
set(gca, 'Projection','perspective')
axis equal
axis manual

figure()
f = Eperp(s,[1;0;0])*r;
f = subs(f, [ty tz], [u v]);
fsurf(f(1),f(2),f(3),[-pi/2 pi/2 -pi pi],'FaceAlpha',0.5);
%f = Eperp(s,[1;0;0]);
%fsurf(@(u,v)u,@(u,v)v,f,[0 2*pi/3 -pi/2 pi/2]);
title('Eperp')
set(gca, 'Projection','perspective')
axis equal
axis manual

figure()
f = Sparr(s,[1;0;0])*r;
f = subs(f, [ty tz], [u v]);
fsurf(f(1),f(2),f(3),[-pi/2 pi/2 -pi pi],'FaceAlpha',0.5);
%f = Sparr(s,[1;0;0]);
%fsurf(@(u,v)u,@(u,v)v,f,[0 2*pi/3 -pi/2 pi/2]);
title('Sparr')
set(gca, 'Projection','perspective')
axis equal
axis manual

figure()
f = Sperp(s,[1;0;0])*r;
f = subs(f, [ty tz], [u v]);
fsurf(f(1),f(2),f(3),[-pi/2 pi/2 -pi pi],'FaceAlpha',0.5);
%f = Sperp(s,[1;0;0]);
%fsurf(@(u,v)u,@(u,v)v,f,[0 2*pi/3 -pi/2 pi/2]);
title('Sperp')
set(gca, 'Projection','perspective')
axis equal
axis manual