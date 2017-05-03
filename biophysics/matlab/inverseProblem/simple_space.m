%function SS = simple_space()

Ndipoles_per_level = [10, 10];
levelDepths = [0.060, 0.040];

rSensors = 0.100; % mm
rSources = 0.060; dSources = sqrt(0.5*rSources^2);

xx = linspace(-dSources,dSources,Ndipoles_per_level(1));
yy = linspace(-dSources,dSources,Ndipoles_per_level(2));

%sourceSpace = zeros(Ndipoles_per_level(1), Ndipoles_per_level(2), length(levelDepths));
sourceLoc = zeros(prod([length(levelDepths), Ndipoles_per_level]), 3);
for kk = 1:length(levelDepths)
    curind = (kk-1)*prod(Ndipoles_per_level) + 1;
    for jj = 1:Ndipoles_per_level(2)
        for ii = 1:Ndipoles_per_level(1)
            sourceLoc(curind,:) = [xx(ii), yy(jj), levelDepths(kk)];
            curind = curind + 1;
        end
    end
end
nSources = size(sourceLoc,1);

clf; hold on
hhss = plot3(sourceLoc(:,1),sourceLoc(:,2),sourceLoc(:,3),'b.');
set(hhss,'MarkerSize',20);


nodefile = 'icosph4.grid';
trifile = 'icosph4.tri';
posBn = rSensors*load(nodefile);
triBn = load(trifile);
[posBn, triBn] = limitTRI(posBn,triBn); 
%trimesh(triBn, posBn(:,1),posBn(:,2),posBn(:,3)); hold on

nodefile = 'icosph3.grid';
trifile = 'icosph3.tri';

posBn = rSensors*load(nodefile);
triBn = load(trifile);

[sensorPos, sensorTri] = limitTRI(posBn,triBn); 
nSensors = size(sensorPos,1);
hh = plot3(sensorPos(:,1),sensorPos(:,2),sensorPos(:,3),'k.');
set(hh,'MarkerSize',20);
axis equal; axis tight; hold on
% %shading interp;%  colorbar
view(2)
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')

%% calculate lead field, aka. gain matrix
G=zeros(size(sensorPos,1),nSources*2);
for ii = 1:nSensors
    tmp = zeros(2,nSources);
    for jj = 1:nSources
        tmp(1,jj) = magsphere_Sarvas(sourceLoc(jj,:), [1,0,0], sensorPos(ii,:)); 
        tmp(2,jj) = magsphere_Sarvas(sourceLoc(jj,:), [0,1,0], sensorPos(ii,:)); 
        G(ii,2*(jj-1)+1) = magsphere_Sarvas(sourceLoc(jj,:), [0,1,0], sensorPos(ii,:)); 
    end
    G(ii,:) = reshape(tmp, 1, nSources*2);
end



% [sensorPos,sensorTri] = load_sensors();
% 
% hh = plot3(sensorPos(:,1),sensorPos(:,2),sensorPos(:,3),'k.');
% axis equal; axis tight;
% set(hh,'MarkerSize',12);
% %shading interp;%  colorbar
% 
% trimesh(sensorTri,sensorPos);
% 
% function [sensorPos, sensorTri] = load_sensors()
% 
% load SensorsM
% 
% nsensors = size(sensLoc, 1);
% nchips = nsensors/3;
% sensorPos = zeros(nchips, 3);
% for ii=1:nchips,
%     
%     sensorPos(ii,:) = sensLoc((ii-1)*3+1,1:3);
%     
% end
% sensorTri = delaunay(sensorPos);

