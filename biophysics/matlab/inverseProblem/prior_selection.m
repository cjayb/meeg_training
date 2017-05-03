%function SS = simple_space()

Ndipoles_per_level = [3, 3];
levelDepths = [0.060];

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

clf; hold on
hhss = plot3(sourceLoc(:,1),sourceLoc(:,2),sourceLoc(:,3),'b.');
set(hhss,'MarkerSize',20);


% nodefile = 'icosph4.grid';
% trifile = 'icosph4.tri';
% posBn = rSensors*load(nodefile);
% triBn = load(trifile);
% [posBn, triBn] = limitTRI(posBn,triBn); 
% %trimesh(triBn, posBn(:,1),posBn(:,2),posBn(:,3)); hold on

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
%fld1 = 2*Gm(:,2) + 2.0*Gm(:,9) + 2*Gm(:,18);
q1   = zeros(size(Gm,2),1); q1(2) = 2;q1(9) = 2; q1(18) = 2;
%fld2 = 0.2*sum(Gm(:,comps),2) + 0.2*sum(Gm(:,comps(evn)),2);

q2   = zeros(size(Gm,2),1);
comps = [1:4, 7:12, 15:18];
evn = ~mod(comps,2);
q2(comps) = 0.275;
q2(comps(evn)) = q2(comps(evn)) + 0.275;
axlim = [-10, 13];
S1 = 0.9; S2 = 0.2;
%%
q1 = zeros(size(Gm,2),1); q1([3,15]) = [1.2,1.2]; q1(9) = - 0.4;
q2 = zeros(size(Gm,2),1); q2([3,15]) = [0.5,0.5]; q2(9) = 1.;
axlim = [-10, 10];
S1 = 1.2; S2 = 0.4;

%%
fld1 = Gm*q1*1e7;%*1e-8*1e15 -> Units will be in fT, assuming q in units of 10 nAm (weird, consider changing)!
fld2 = Gm*q2*1e7;

%%

figure(101);clf; 
subplot(2,2,1);
plot_fieldpattern3D(sensorPos,sensorTri,fld1,sensorPos); cb(1) = colorbar;
subplot(2,2,2);
plot_fieldpattern3D(sensorPos,sensorTri,fld2,sensorPos); cb(2) = colorbar;
subplot(2,2,3);
plot_sourcepattern2D(sourceLoc, q1,S1); axis equal; set(gca,'Xlim', axlim*1e-3); set(gca,'Ylim', axlim*1e-3)
subplot(2,2,4);
plot_sourcepattern2D(sourceLoc, q2,S2); axis equal; set(gca,'Xlim', axlim*1e-3); set(gca,'Ylim', axlim*1e-3);
colormap(gray(256))


n1_L2 = norm(q1)*10;
n2_L2 = norm(q2)*10;
n1_L1 = norm(q1,1)*10;
n2_L1 = norm(q2,1)*10;

fprintf(1,'\tL2\tL1\n')
fprintf(1,'Q1\t%.2f\t%.2f\n',n1_L2,n1_L1);
fprintf(1,'Q2\t%.2f\t%.2f\n',n2_L2,n2_L1);


