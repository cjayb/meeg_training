clear

boundary_cols = [0.5,0.5,0;...
    0,0.5,0.5];


R = [8.0, 8.1, 8.6, 9.2] / 100; %m
sigmas = [0.33, 1, 0.0042, 0.33];
center = [0,0,0]; % MUST be origin!

thetaEl = [-pi/8; pi/8];
Re = R(end) * [zeros(size(thetaEl)), sin(thetaEl), cos(thetaEl)];

% how many grid points
N = 128;
ylim = [-1, 1];
zlim = [-1, 1];
gpy = round(diff(ylim)/2)*N; % even to avoid origin in set
gpz = round(diff(zlim)/2)*N; % even to avoid origin in set
yy=linspace(ylim(1),ylim(2),gpy) * R(1);   
zz=linspace(zlim(1),zlim(2),gpz) * R(1);   

[Y,Z]=meshgrid(yy,zz); 
% Y = reshape(Y, N^2,1);
% Z = reshape(Z, N^2,1);
r = sqrt(Y.^2 + Z.^2);
inInd = r <= 1.05*R(1); % a bit larger than necc....
inInd_alpha = r <= R(1); % tight

% Y = Y(r < R(1));
% Z = Z(r < R(1));

Rq = [zeros(size(Y(inInd(:)))), Y(inInd(:)), Z(inInd(:))];

% figure(1); clf;
% plot(Y(inInd(:)),Z(inInd(:)),'.');hold on
% plot(Re(:,2), Re(:,3), 'ro')

G = bst_eeg_sph_cjb(Rq, Re, center, R, sigmas);

EEG = reshape(G(1,:), 3, sum(inInd(:))) - reshape(G(2,:), 3, sum(inInd(:)));

C = zeros(size(Y)); C(inInd(:)) = sqrt(dot(EEG,EEG));
figure(2); clf
%scatter(Y,Z,5,sqrt(dot(EEG,EEG)))

C = C/max(C(inInd_alpha(:)));

steplist = [0, 0.1, 0.25, 0.5, 0.75, 0.9];
contsteps = steplist;
% contsteps = max(C(inInd_alpha(:))) * steplist;

load parula_cmap.mat
cmap(end,:) = [1,1,1];
circleMask = contsteps(end)*ones(size(r));
circleAlpha = ones(size(r));
circleAlpha(inInd_alpha) = 0;

[cc,ch] = contourf(Y,Z,C, contsteps, 'LineStyle','-.', 'HitTest','off'); 
colormap(cmap); colorbar

% lh=clabel(cc, [0.1,0.5], 'FontSize',14,'Color','w');

set(ch,'TextList', [0.5])
set(ch,'LineColor', 'w')
set(ch,'ShowText', 'on')

axis xy; axis equal; hold on
imagesc(yy,zz,circleMask,...
    'AlphaData', circleAlpha, 'HitTest','off');

nTheta = 100;
theta = linspace(0,2*pi,nTheta);
plot(R(1)*cos(theta), R(1)*sin(theta), 'Color',boundary_cols(1,:), ...
    'LineWidth', 2, 'Hittest', 'off');
plot(R(end)*cos(theta), R(end)*sin(theta), 'Color',boundary_cols(2,:), ...
    'LineWidth', 4, 'Hittest', 'off');
plot(Re(:,2), Re(:,3), 'ro', 'MarkerSize',10, 'MarkerFaceColor', 'r')

%% 3D EEG plot
w = [-1, 1];
w = w/length(w);
er = zeros(size(w));
senstype = 'EEG';
figure(21)
plot_lead_on_3Dsph(Re, w, er, R, sigmas, senstype, cmap)

%% old-school magsphere_Sarvas.m

% Rm = Re;
% %Rq = [0, 0.05,0.05; 0, -0.05,0.05];
% 
% % only the x-component contributes in yz-plane
% Gm = zeros(size(Re,1), size(Rq,1));
% q = [1,0,0]; %q = eye(3);
% totcalcs = prod(size(Gm));
% 
% for ii = 1:size(Rm, 1)
%     for jj = 1:size(Rq, 1)
%         if ~mod(jj,1000)
%             fprintf(1, '%d of %d\n',(ii-1)*size(Rq, 1) + jj, totcalcs)
%         end
% %         for kk = 1:3
% %             fprintf(1, 'ii=%d, Gm2=%d:%d\n', ii,...
% %                 9*(jj-1) + (kk-1)*3+1, 9*(jj-1) + (kk-1)*3+3);
% %             Gm(ii,9*(jj-1) + (kk-1)*3+1:9*(jj-1) + (kk-1)*3+3) = ...
% %                 magsphere_Sarvas(Rq(jj,:), q(kk,:), Rm(ii,:));
%         Gm(ii,jj) = ...
%             magsphere_Sarvas(Rq(jj,:), q, Rm(ii,:));
% %         end
%     end
% end

%% vectorised!
%thetaEl = thetaEl(1);
Rm = (R(end)+0.03) * [zeros(size(thetaEl)), sin(thetaEl), cos(thetaEl)];

% Rm = Re;
% only the x-component contributes in yz-plane
% for all 3 q's, loop over them
Gm = zeros(size(Rm,1), size(Rq,1));
q = [1,0,0]; %q = eye(3);

for ii = 1:size(Rm, 1)
        Gm(ii,:) = ...
            magsphere_Sarvas_vec(Rq, q, Rm(ii,:));
end

%%
% MEG = reshape(Gm(1,:), 3, sum(inInd(:))) + reshape(Gm(2,:), 3, sum(inInd(:)));
% MEG = abs(Gm(1,:)) + abs(Gm(2,:));
MEG = sqrt(sum(Gm.^2,1));
Cm = zeros(size(Y)); Cm(inInd(:)) = MEG;

Cm = Cm/max(Cm(inInd_alpha(:)));
steplist = [0, 0.1, 0.25, 0.5, 0.75, 0.9];

figure(3);clf; 
contsteps = steplist;
% contsteps = max(Cm(inInd_alpha(:))) * [0:0.1:0.5, 0.8, 0.99];

load parula_cmap.mat
cmap(end,:) = [1,1,1];
circleMask = contsteps(end)*ones(size(r));
circleAlpha = ones(size(r));
circleAlpha(inInd_alpha) = 0;

[cc,ch] = contourf(Y,Z,Cm, contsteps, 'LineStyle','-.', 'HitTest','off'); 
colormap(cmap); colorbar
axis xy; axis equal; hold on
imagesc(yy,zz,circleMask,...
    'AlphaData', circleAlpha, 'HitTest','off');

nTheta = 100;
theta = linspace(0,2*pi,nTheta);
plot(R(1)*cos(theta), R(1)*sin(theta), 'Color',boundary_cols(1,:), ...
    'LineWidth', 2, 'Hittest', 'off');
plot(R(end)*cos(theta), R(end)*sin(theta), 'Color',boundary_cols(2,:), ...
    'LineWidth', 4, 'Hittest', 'off');
plot(Rm(:,2), Rm(:,3), 'bo', 'MarkerSize',10, 'MarkerFaceColor', 'c')

% lh=clabel(cc, [0.1,0.5], 'FontSize',14,'Color','r');
set(ch,'TextList', [0.5])
set(ch,'LineColor', 'w')
set(ch,'ShowText', 'on')

%% now for some real sensors: mags
% deal with one sensor at a time
% thetaMagPos = thetaEl(1); % zero-angle is on z-axis, opens clockwise
thetaMagPos = 0; % zero-angle is on z-axis, opens clockwise
thetaMagRot = -thetaMagPos; % zero-angle is on y-axis, opens counter-cl.
senstype = 'ygrad';
% senstype = 'mag';
[rC, w, er] = coildefs(senstype);

% first rotate in place
rotMat = [1,            0,             0; ...
          0, cos(thetaMagRot), -sin(thetaMagRot); ...
          0, sin(thetaMagRot), cos(thetaMagRot)];

rC = (rotMat * rC')';
er = (rotMat * er')';

% move center to coil position
rHelmet = R(end)+0.03;
Rm = rHelmet * [zeros(size(thetaMagPos)), sin(thetaMagPos), cos(thetaMagPos)];

rC = rC + repmat(Rm, size(rC, 1), 1);
% plot(rC(:,2), rC(:,3),'ro')

coilframe = coil_outline(Rm, rotMat, senstype);

%%
Gmag = zeros(3, size(rC,1), size(Rq,1));
q = eye(3); % now there's an x-component to the r-vectors

for ii = 1:size(rC, 1)
    for kk = 1:3
        % bring in the weights here for easy summing later
        Gmag(kk,ii,:) = w(ii) * ...
            magsphere_Sarvas_vec(Rq, q(kk,:), rC(ii,:), er(ii,:));
    end
end
% take the norm over q-orientations and sum (weigths applied above)!
% Gmag = squeeze(sqrt(dot(Gmag,Gmag,1)));
% deal with w's
% MEG = sum(Gmag.*repmat(w',1, size(Gmag,2)), 1);

% sum first, then norm!
Gmag = squeeze(sum(Gmag, 2)); % sum over integration points
MEG = squeeze(sqrt(dot(Gmag,Gmag,1)));

% MEG = reshape(Gm(1,:), 3, sum(inInd(:))) + reshape(Gm(2,:), 3, sum(inInd(:)));
%EEG = reshape(G(1,:), 3, sum(inInd(:))) - reshape(G(2,:), 3, sum(inInd(:)));
%MEG = mean(Gmag, 1);

Cm = zeros(size(Y)); Cm(inInd(:)) = MEG;

Cm = Cm/max(Cm(inInd_alpha(:)));
steplist = [0, 0.1, 0.25, 0.5, 0.75, 0.9];

figure(4);clf; 
contsteps = steplist;
% contsteps = max(Cm(inInd_alpha(:))) * [0:0.1:0.5, 0.8, 0.99];

load parula_cmap.mat
cmap(end,:) = [1,1,1];
circleMask = contsteps(end)*ones(size(r));
circleAlpha = ones(size(r));
circleAlpha(inInd_alpha) = 0;

[cc,ch] = contourf(Y,Z,Cm, contsteps, 'LineStyle','-.', 'HitTest','off'); 
colormap(cmap); colorbar
axis xy; axis equal; hold on
imagesc(yy,zz,circleMask,...
    'AlphaData', circleAlpha, 'HitTest','off');

nTheta = 100;
theta = linspace(0,2*pi,nTheta);
plot(R(1)*cos(theta), R(1)*sin(theta), 'Color',boundary_cols(1,:), ...
    'LineWidth', 2, 'Hittest', 'off');
plot(R(end)*cos(theta), R(end)*sin(theta), 'Color',boundary_cols(2,:), ...
    'LineWidth', 4, 'Hittest', 'off');
plot(Rm(:,2), Rm(:,3), 'bo', 'MarkerSize',10, 'MarkerFaceColor', 'c')

% lh=clabel(cc, [0.1,0.5], 'FontSize',14,'Color','r');
set(ch,'TextList', [0.5])
set(ch,'LineColor', 'w')
set(ch,'ShowText', 'on')

figure(41); clf
plot_lead_on_3Dsph(rC, w, er, R, sigmas, senstype, cmap)
hold on
sh=plot3(coilframe(:,1),coilframe(:,2),coilframe(:,3));
% coilframe = coil_outline(Rm, rotMat, 'outline');
% sh=plot3(coilframe(:,1),coilframe(:,2),coilframe(:,3),'r');

% sh=plot3(coilframe(1,1),coilframe(1,2),coilframe(1,3),'ro');
% sh=plot3(coilframe(end,1),coilframe(end,2),coilframe(end,3),'go');
% sh=plot3(coilframe(8,1),coilframe(8,2),coilframe(8,3),'cx');
% sh=plot3(coilframe(9,1),coilframe(9,2),coilframe(9,3),'cx');
% sh=plot3(coilframe(10,1),coilframe(10,2),coilframe(10,3),'cx');
% sh=plot3(coilframe(end-3:end,1),coilframe(end-3:end,2),coilframe(end-3:end,3),'cx');
% set(sh,'FaceColor','r')
% set(sh,'FaceAlpha',0.3)
% set(sh,'EdgeColor','r')
%% now for some real sensors: mags
% deal with one sensor at a time
figure(5);clf; 
sensors = {'mag','xgrad','ygrad'};

for ss = 1:3
    subplot(1,3,ss)
    
%     thetaMagPos = thetaEl(1); % zero-angle is on z-axis, opens clockwise
    thetaMagPos = 0; % zero-angle is on z-axis, opens clockwise
    thetaMagRot = -thetaMagPos; % zero-angle is on y-axis, opens counter-cl.
    [rC, w, er] = coildefs(sensors{ss});
    
    % first rotate in place
    rotMat = [1,            0,             0; ...
        0, cos(thetaMagRot), -sin(thetaMagRot); ...
        0, sin(thetaMagRot), cos(thetaMagRot)];
    
    rC = (rotMat * rC')';
    er = (rotMat * er')';
    
    % move center to coil position
    rHelmet = R(end)+0.03;
    Rm = rHelmet * [zeros(size(thetaMagPos)), sin(thetaMagPos), cos(thetaMagPos)];
    
    rC = rC + repmat(Rm, size(rC, 1), 1);
    plot(rC(:,2), rC(:,3),'ro')
    
    %%
    Gmag = zeros(3, size(rC,1), size(Rq,1));
    q = eye(3); % now there's an x-component to the r-vectors
    
    for ii = 1:size(rC, 1)
        for kk = 1:3
            % bring in the weights here for easy summing later
            Gmag(kk,ii,:) = w(ii) * ...
                magsphere_Sarvas_vec(Rq, q(kk,:), rC(ii,:), er(ii,:));
        end
    end
    % take the norm over q-orientations and sum (weigths applied above)!
    % Gmag = squeeze(sqrt(dot(Gmag,Gmag,1)));
    % deal with w's
    % MEG = sum(Gmag.*repmat(w',1, size(Gmag,2)), 1);
    
    % sum first, then norm!
    Gmag = squeeze(sum(Gmag, 2)); % sum over integration points
    MEG = squeeze(sqrt(dot(Gmag,Gmag,1)));
    %%
    % MEG = reshape(Gm(1,:), 3, sum(inInd(:))) + reshape(Gm(2,:), 3, sum(inInd(:)));
    %EEG = reshape(G(1,:), 3, sum(inInd(:))) - reshape(G(2,:), 3, sum(inInd(:)));
    %MEG = mean(Gmag, 1);
    
    Cm = zeros(size(Y)); Cm(inInd(:)) = MEG;
    
    Cm = Cm/max(Cm(inInd_alpha(:)));
    steplist = [0, 0.1, 0.25, 0.5, 0.75, 0.9];
    
    contsteps = steplist;
    % contsteps = max(Cm(inInd_alpha(:))) * [0:0.1:0.5, 0.8, 0.99];
    
    load parula_cmap.mat
    cmap(end,:) = [1,1,1];
    circleMask = contsteps(end)*ones(size(r));
    circleAlpha = ones(size(r));
    circleAlpha(inInd_alpha) = 0;
    
    [cc,ch] = contourf(Y,Z,Cm, contsteps, 'LineStyle','-.', 'HitTest','off');
    colormap(cmap); colorbar
    axis xy; axis equal; hold on
    imagesc(yy,zz,circleMask,...
        'AlphaData', circleAlpha, 'HitTest','off');
    
    nTheta = 100;
    theta = linspace(0,2*pi,nTheta);
    plot(R(1)*cos(theta), R(1)*sin(theta), 'Color',boundary_cols(1,:), ...
        'LineWidth', 2, 'Hittest', 'off');
    plot(R(end)*cos(theta), R(end)*sin(theta), 'Color',boundary_cols(2,:), ...
        'LineWidth', 4, 'Hittest', 'off');
    plot(Rm(:,2), Rm(:,3), 'bo', 'MarkerSize',10, 'MarkerFaceColor', 'c')
    
    % lh=clabel(cc, [0.1,0.5], 'FontSize',14,'Color','r');
    set(ch,'TextList', [0.5])
    set(ch,'LineColor', 'w')
    set(ch,'ShowText', 'on')
end