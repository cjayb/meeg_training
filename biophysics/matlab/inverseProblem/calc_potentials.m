function V = calc_potentials(sourcePos, sourceMom, elPos, bLead)
% sourcePos (m): Nx3
% sourceMom (A m): Nx3
% elPos (m): Mx3
% bLead (bool): return lead field matrix
%
% V (V): Mx1
if nargin < 4
    bLead = 0;
end

nSensors = size(elPos,1);
nSources= size(sourcePos,1);
V = zeros(nSensors,nSources);
for ii = 1:nSensors
    for jj = 1:nSources
        r = 100*sourcePos(jj,:);  % cm
        m = 1e5*sourceMom(jj,:);  % ->mA cm = x1000 X100
        s = 100*elPos(ii,:);  % cm
        V(ii, jj) = 1e-3*fastpotential_sun(r, m, s);  % mV -> V
    end
end
if ~bLead,
    V = sum(V, 2);
end