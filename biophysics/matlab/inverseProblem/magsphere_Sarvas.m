function Br = magsphere_Sarvas(r0, q, r)

a = r-r0;


R = sqrt(dot(r,r));
A = sqrt(dot(a,a));

F = A*(R*A + R*R - dot(r0,r));


fact1 = (A*A/R + dot(a,r)/A + 2.0*A + 2.0*R);
fact2 = (A + 2.0*R + dot(a,r)/A);


nablaF = fact1*r - fact2*r0;

Qr0 = cross(q, r0);
Qr0r = dot(Qr0, r);

%/*MU4PIF = PERM/(4*PI*F*F);*/
%/*mu = 4*pi*10e-7 => mu/4pi = 1e-7 = 1/10000000*/
MU4PIFinv = F*F*10000000.0;

B = (F*Qr0 - Qr0r*nablaF)/MU4PIFinv;

er = r./norm(r);
Br = dot(B,er);