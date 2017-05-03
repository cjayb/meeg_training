function Br=fastmeg(dip,mom)

nodefile = 'icosph4.grid';
trifile = 'icosph4.tri';
elec = load(nodefile);
tri = load(trifile);

headr = 0.080;  %m
elec = headr*elec; 

Br=zeros(size(elec,1),1);
for ii = 1:size(elec,1)
    Br(ii) = magsphere_Sarvas(dip, mom, elec(ii,:)); % -> V
end

elec = elec*1000;   % m -> mm
trisurf(tri, elec(:,1),elec(:,2),elec(:,3),Br);
axis equal; axis tight;
shading interp;%  colorbar
        
view(2)
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
