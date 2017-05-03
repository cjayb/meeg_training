function V = fasteeg(dip, mom)
%
% USAGE: V = fasteeg(elec, dip, mom, tri),
%   V:      potential (V)
%   elec:   electrode {x,y,z}-coordinates (m)
%   dip:    dipole {x,y,z}-coordinates (m)
%   mom:    dipole {x,y,z}-moments (A m) 
%   tri:    matrix forming triangles for electrodes (OPTIONAL!)
%
%/**********************************************************************/
%/* Sun M., "An Efficient Algorithm for Computing Multishell Spherical */
%/* Volume Conductor Models in EEG Dipole Source Localization",        */
%/* IEEE Trans. Biomed. Eng., 44(12), pp. 1243--52, 1997.              */
%/*                                                                    */
%/* Code copied directly from page 1248 of above publication!          */
%/*                                                                    */
%/* NB! Code expects and returns the following units:                  */
%/*     Dipole location    (r): cm                                     */
%/*     Dipole moment      (m): mA cm                                  */
%/*     Electrode location (s): cm                                     */
%/*     *** Output voltage    : mV                                     */
%/*                                                                    */
%/* The (default) 4-shell conductor model is by Stok (86):             */
%/*                                                                    */
%/*               Conductivity (S/m)     Radius (relative to scalp, R) */
%/*     Brain          0.33                0.8400  (63 mm)             */
%/*     CSF            1.0                 0.8667  (65 mm)             */
%/*     Skull          0.0042              0.9467  (71 mm)             */
%/*     Scalp          0.33                1.0000  (75 mm)             */
%/*                                                                    */
%/*                                                                    */
%/* C. Bailey, 17/06/02                                                */
%/**********************************************************************/

nodefile = 'icosph4.grid';
trifile = 'icosph4.tri';
elec = load(nodefile);
tri = load(trifile);

headr = 0.080;  %m
elec = headr*elec; 

elec = elec*100;   % m -> cm
dip  = dip*100;    % m -> cm
mom  = mom*1e5;    % Am -> mA cm

V=zeros(size(elec,1),1);
for ii = 1:size(elec,1)
    V(ii) = fastpotential_sun(dip, mom, elec(ii,:))/1000; % -> V
end

elec = elec*10;   % cm -> mm
trisurf(tri, elec(:,1),elec(:,2),elec(:,3),V);
axis equal; axis tight;
shading interp;% colorbar

view(2)    
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
