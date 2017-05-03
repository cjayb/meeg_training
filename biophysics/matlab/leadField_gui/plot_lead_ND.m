function plot_lead_ND(handles, rSens, wSens, erSens, R, sigmas, senstype, b2Dquiver)

% assume cm!
rSens = rSens/100;
R = R/100;

steplist = [0, 0.1, 0.25, 0.5, 0.75, 0.9];
contsteps = steplist;
center = [0,0,0];
% colormap
% load parula_cmap.mat
cmap = viridis();
cmap(end,:) = [1,1,1];
quiver_col = [0.0, 0.75, 0.85];

boundary_cols = [0.5,0.5,0;...
    0,0.5,0.5];


% how many grid points
N = 128;
ylim = [-1, 1];
zlim = [-1, 1];
gpy = round(diff(ylim)/2)*N; % even to avoid origin in set
gpz = round(diff(zlim)/2)*N; % even to avoid origin in set
yy=linspace(ylim(1),ylim(2),gpy) * R(1);   
zz=linspace(zlim(1),zlim(2),gpz) * R(1);   

[Y,Z]=meshgrid(yy,zz); 
r = sqrt(Y.^2 + Z.^2);
inInd = r <= 1.05*R(1); % a bit larger than necc....
inInd_alpha = r <= R(1); % tight

Rq = [zeros(size(Y(inInd(:)))), Y(inInd(:)), Z(inInd(:))];

if strcmp(senstype(1:3), 'EEG')
    
    G = bst_eeg_sph_cjb(Rq, rSens, center, R, sigmas);
    
%     REF = reshape(G(1,:), 3, sum(inInd(:)));
    EEG = zeros(3, size(rSens,1), size(Rq,1));
    for ii = 1:size(rSens,1)
%     for ii = 2:size(rSens,1)
%         EEG(:,ii-1,:) = reshape(G(ii,:), 3, sum(inInd(:))) - REF;
        EEG(:,ii,:) = wSens(ii) * reshape(G(ii,:), 3, sum(inInd(:)));
    end
    
%     EEG = squeeze(sum(EEG,2));
%     EEG = sqrt(dot(EEG,EEG, 1));

    if size(rSens, 1) > 2
        b2Dquiver = 0;
    else
        Gvec = squeeze(sum(EEG, 2)); % sum over coils, for quiver
    end
    EEG = squeeze(sqrt(dot(EEG,EEG, 1)));
    EEG = squeeze(sum(EEG,1));
    C = zeros(size(Y)); C(inInd(:)) = EEG;


else
    nSensors = size(rSens,1)/4;
    G = zeros(3, nSensors, 4, size(Rq,1));
    q = eye(3); % now there's an x-component to the r-vectors
    
%     for ii = 1:size(rSens, 1)
    for ii = 1:nSensors
        for jj = 1:4
            for kk = 1:3
                % bring in the weights here for easy summing later
                G(kk,ii,jj,:) = wSens((ii-1)*4 + jj) * ...
                    magsphere_Sarvas_vec(Rq, q(kk,:), ...
                    rSens((ii-1)*4 + jj,:), erSens((ii-1)*4 + jj,:));
                %                 magsphere_Sarvas_vec(Rq, q(kk,:), rSens(ii,:), erSens(ii,:));
            end
        end
    end
%     % sum first, then norm!
%     G = squeeze(sum(G, 2)); % sum over integration points
%     MEG = squeeze(sqrt(dot(G,G,1)));
% %     MEG = squeeze(sqrt(dot(G,G,1)));
% %     MEG = squeeze(sum(MEG, 1)); % sum over integration points
% %     %
% %     MEG = squeeze(sum(MEG, 1));

    G = squeeze(sum(G, 3)); % sum over integration points
    if nSensors == 1,
        Gvec = G;
    else
        b2Dquiver = 0; % no sense to show sum vector...
%         Gvec = squeeze(sum(G, 2)); % sum over coils, for quiver
    end
    
    MEG = squeeze(sqrt(dot(G,G,1))); % norm over q
    MEG = squeeze(sum(MEG, 1)); % sum over coils
    
    C = zeros(size(Y)); C(inInd(:)) = MEG;
    
   
end

C = C/max(C(inInd_alpha(:)));

axes(handles.axes1)
children = get(handles.axes1,'Children');
for jj = 1:length(children)
    % there's no contour-type in 2014a! add test for hggroup!
    if strcmp(get(children(jj), 'Type'), 'contour') || ...
            strcmp(get(children(jj), 'Type'), 'quiver') || ...
            strcmp(get(children(jj), 'Type'), 'hggroup') || ... 
            strcmp(get(children(jj), 'Type'), 'image')
        delete(children(jj))
    end
end


[cc,ch] = contourf(Y*100,Z*100,C, contsteps, 'LineStyle','-.', 'HitTest','off');
colormap(cmap); colorbar
axis xy; axis equal; hold on
xlabel('posterior-anterior'); ylabel('inferior-superior');

set(ch,'TextList', [0.5])
set(ch,'LineColor', 'w')
set(ch,'ShowText', 'on')

if b2Dquiver
    skipafew = 45;
    red_Gvec = Gvec(:, 1:skipafew:end);
    len_Gvec = sqrt(sum(red_Gvec.*red_Gvec, 1));
    normGvec = red_Gvec ./ repmat(len_Gvec, 3, 1);
    qh=quiver(100*Rq(1:skipafew:end,2),100*Rq(1:skipafew:end,3),...
        normGvec(2,:)',normGvec(3,:)');
    set(qh,'Color',quiver_col, 'LineWidth', 0.5, 'AutoScale','on')
end

circleMask = contsteps(end)*ones(size(r));
circleAlpha = ones(size(r));
circleAlpha(inInd_alpha) = 0;
imagesc(yy*100,zz*100,circleMask,...
    'AlphaData', circleAlpha, 'HitTest','off');

children = get(handles.axes1,'Children');
if b2Dquiver
    children = [children(4:end); children(1:3)];
else
    children = [children(3:end); children(1:2)];
end
set(handles.axes1,'Children', children);
% disp(children)

% add_circle(R(1), '-', 2, boundary_cols(1,:))
% add_circle(R(end), '-', 4, boundary_cols(2,:))
% nTheta = 100;
% theta = linspace(0,2*pi,nTheta);
% plot(R(1)*cos(theta), R(1)*sin(theta), 'Color',boundary_cols(1,:), ...
%     'LineWidth', 2, 'Hittest', 'off');
% plot(R(end)*cos(theta), R(end)*sin(theta), 'Color',boundary_cols(2,:), ...
%     'LineWidth', 4, 'Hittest', 'off');

%% 3D

axes(handles.axes2); cla

nodefile = 'icosph4.grid';
trifile = 'icosph4.tri';
Rq = load(nodefile)*R(1); %innermost sphere
tri = load(trifile);

q = eye(3); 

if strcmp(senstype(1:3), 'EEG')
    %pass
    % here wSens can define Ref!
%     R = [8.0, 8.1, 8.6, 9.2] / 100; %m
%     sigmas = [0.33, 1, 0.0042, 0.33];
    G = zeros(3, size(rSens,1), size(Rq,1));
    
    center = [0,0,0]; % MUST be origin!
    Gflat = bst_eeg_sph_cjb(Rq, rSens, center, R, sigmas);
    for ii = 1:size(Gflat,1)
        G(:, ii, :) = wSens(ii) * reshape(Gflat(ii,:), 3, size(Rq,1));
    end
    faceAlpha = 1;
else
    nSensors = size(rSens,1)/4;
    G = zeros(3, nSensors, 4, size(Rq,1));
    for ii = 1:nSensors
        for jj = 1:4
            for kk = 1:3
                % bring in the weights here for easy summing later
                G(kk,ii,jj,:) = wSens((ii-1)*4 + jj) * ...
                    magsphere_Sarvas_vec(Rq, q(kk,:), ...
                    rSens((ii-1)*4 + jj,:), erSens((ii-1)*4 + jj,:));
            end
        end
    end
    G = squeeze(sum(G, 3)); % sum over integration points

    faceAlpha = 1;
end

% norm, then sum
LF = squeeze(sqrt(dot(G,G,1)));

if numel(G) == 3*size(Rq,1) % only one coil!
    Gvec = G; % sum over sensors
else
    Gvec = squeeze(sum(G, 2)); % sum over sensors
    LF = squeeze(sum(LF,1));
end
% G = squeeze(sum(G, 2)); % sum over sensors
% LF = squeeze(sqrt(dot(G,G,1)));

LF = LF/max(LF(:));

Rq = Rq * 100;

ht = trisurf(tri, Rq(:,1),Rq(:,2),Rq(:,3),LF);
set(ht,'FaceAlpha',faceAlpha)

if b2Dquiver
    hold on
    qh=quiver3(Rq(:,1),Rq(:,2),Rq(:,3),Gvec(1,:)',Gvec(2,:)',Gvec(3,:)');
    hold off
    set(qh,'Color',quiver_col)
end
axis equal; axis tight;
shading interp;
colormap(cmap);  colorbar

view(90,90)
xlabel('left-right')
ylabel('posterior-anterior')
zlabel('inferior-superior')




