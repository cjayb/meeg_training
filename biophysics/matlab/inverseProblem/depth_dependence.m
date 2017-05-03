%function SS = simple_space()

Ndipoles_per_level = [3, 3];
levelDepths = [0.060, 0.045];

rSensors = 0.100; % mm
rSources = 0.010; dSources = sqrt(0.5*rSources^2);

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

% clf; hold on
% hhss = plot3(sourceLoc(:,1),sourceLoc(:,2),sourceLoc(:,3),'b.');
% set(hhss,'MarkerSize',20);

nodefile = 'icosph3.grid';
trifile = 'icosph3.tri';

posBn = rSensors*load(nodefile);
triBn = load(trifile);

[sensorPos, sensorTri] = limitTRI(posBn,triBn); 
nSensors = size(sensorPos,1);
% hh = plot3(sensorPos(:,1),sensorPos(:,2),sensorPos(:,3),'k.');
% set(hh,'MarkerSize',20);
% axis equal; axis tight; hold on
% % %shading interp;%  colorbar
% view(2)
% xlabel('x-axis')
% ylabel('y-axis')
% zlabel('z-axis')

%% calculate lead field, aka. gain matrix
Gm=zeros(size(sensorPos,1),nSources*2);
for ii = 1:nSensors
    tmp = zeros(2,nSources);
    for jj = 1:nSources
        tmp(1,jj) = magsphere_Sarvas(sourceLoc(jj,:), [1,0,0], sensorPos(ii,:)); 
        tmp(2,jj) = magsphere_Sarvas(sourceLoc(jj,:), [0,1,0], sensorPos(ii,:)); 
        %G(ii,2*(jj-1)+1) = magsphere_Sarvas(sourceLoc(jj,:), [0,1,0], sensorPos(ii,:)); 
    end
    Gm(ii,:) = reshape(tmp, 1, nSources*2);
end
%%
q1 = zeros(size(Gm,2),1); q1(9) = 10;
q2 = zeros(size(Gm,2),1); q2(18+9) = 10;
axlim = [-10, 10];
S1 = 1.; S2 = 1.;

%%
fld1 = Gm*q1*1e7;%*1e-8*1e15 -> Units will be in fT!
fld2 = Gm*q2*1e7;

%%

figure(201);clf; 
subplot(2,2,1);
plot_fieldpattern3D(sensorPos,sensorTri,fld1,sensorPos); cb(1) = colorbar;
subplot(2,2,2);
plot_fieldpattern3D(sensorPos,sensorTri,fld2,sensorPos); cb(2) = colorbar;
subplot(2,2,3);
plot_sourcepattern2D(sourceLoc, q1,S1); axis equal; set(gca,'Xlim', axlim*1e-3); set(gca,'Ylim', axlim*1e-3)
subplot(2,2,4);
plot_sourcepattern2D(sourceLoc, q2,S2); axis equal; set(gca,'Xlim', axlim*1e-3); set(gca,'Ylim', axlim*1e-3);
%colormap(gray(256))
%colormap(cmap_freaky)
% colormap(cmap_pastel)
colormap(parula)

n1_L2 = norm(q1);
n2_L2 = norm(q2);
n1_L1 = norm(q1,1);
n2_L1 = norm(q2,1);

fprintf(1,'\tL2\tL1\n')
fprintf(1,'Q1\t%.2f\t%.2f\n',n1_L2,n1_L1);
fprintf(1,'Q2\t%.2f\t%.2f\n',n2_L2,n2_L1);


