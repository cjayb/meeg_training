function Br = magsphere_Sarvas_vec(r0, q, r, er)
% r0: positions of q [N x 3]
%  q: single-orientation q [1 x 3]
%  r: single sensor position [1 x 3]

rN = repmat(r, size(r0,1), 1);
qN = repmat(q, size(r0,1), 1);

a = rN-r0;


R = sqrt(dot(r,r, 2)); % dot in 2nd dim -> [N x 1]
A = sqrt(dot(a,a, 2)); % dot in 2nd dim

F = A.*(R.*A + R.*R - dot(r0,rN, 2));


fact1 = (A.*A./R + dot(a,rN, 2)./A + 2.0*A + 2.0*R);
fact2 = (A + 2.0*R + dot(a,rN, 2)./A);


nablaF = repmat(fact1,1,3).*rN - repmat(fact2,1,3).*r0;

Qr0 = cross(qN, r0, 2);
Qr0r = dot(Qr0, rN, 2);

%/*MU4PIF = PERM/(4*PI*F*F);*/
%/*mu = 4*pi*10e-7 => mu/4pi = 1e-7 = 1/10000000*/
MU4PIFinv = F.*F*10000000.0;

B = (repmat(F,1,3).*Qr0 - repmat(Qr0r,1,3).*nablaF)./...
    repmat(MU4PIFinv,1,3);

% assume point-magnetometer at r
if nargin < 4
    er = r./norm(r);
end
ER = repmat(er,size(r0,1), 1); 
Br = dot(B,ER,2);