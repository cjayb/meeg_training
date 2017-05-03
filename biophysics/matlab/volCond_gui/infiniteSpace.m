function V = infiniteSpace(rq, q, R, M, Y, Z)

Rq = sqrt(dot(rq,rq));
Q = sqrt(dot(q,q));

lam1 = 0; % do this in 2D, let qx==0
mu1 = rq(2)/Rq;
nu1 = rq(3)/Rq; 

psix = q(1)/Q;
psiy = q(2)/Q;
psiz = q(3)/Q; 

fR = Rq/R;

r = sqrt(Y.^2 + Z.^2);
lam = zeros(size(r));
mu = cos(atan2(Z,Y)); % 
%nu = cos(atan2(Z,Y) - pi/2); % WHY MINUS ?! PI/2?!
nu = cos(atan2(-Y,Z)); % 

gamma = lam1*lam + mu1*mu + nu1*nu;
rpq = sqrt(r.^2 + (fR)^2 - 2*fR*r.*gamma);

V = M*( psix*(r.*lam -fR*lam1) + psiy*(r.*mu -fR*mu1) + ...
    psiz*(r.*nu -fR*nu1) ) ./ rpq.^3;

