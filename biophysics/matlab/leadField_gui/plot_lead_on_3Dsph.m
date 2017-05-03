function plot_lead_on_3Dsph(rSens, wSens, erSens, R, sigmas, senstype, cmap)

nodefile = 'icosph4.grid';
trifile = 'icosph4.tri';
Rq = load(nodefile)*R(1); %innermost sphere
tri = load(trifile);

G = zeros(3, size(rSens,1), size(Rq,1));
q = eye(3); 

if strcmp(senstype, 'EEG')
    %pass
    % here wSens can define Ref!
%     R = [8.0, 8.1, 8.6, 9.2] / 100; %m
%     sigmas = [0.33, 1, 0.0042, 0.33];
    center = [0,0,0]; % MUST be origin!
    Gflat = bst_eeg_sph_cjb(Rq, rSens, center, R, sigmas);
    for ii = 1:size(Gflat,1)
        G(:, ii, :) = wSens(ii) * reshape(Gflat(ii,:), 3, size(Rq,1));
    end
    faceAlpha = 0.5;
else
    for ii = 1:size(rSens, 1)
        for kk = 1:3
            % bring in the weights here for easy summing later
            G(kk,ii,:) = wSens(ii) * ...
                magsphere_Sarvas_vec(Rq, q(kk,:), rSens(ii,:), erSens(ii,:));
        end
    end
    faceAlpha = 1;
end
% take the norm over q-orientations and sum (weigths applied above)!
% Gmag = squeeze(sqrt(dot(Gmag,Gmag,1)));
% deal with wSens's
% MEG = sum(Gmag.*repmat(wSens',1, size(Gmag,2)), 1);

% sum first, then norm!
G = squeeze(sum(G, 2)); % sum over integration points
LF = squeeze(sqrt(dot(G,G,1)));

ht = trisurf(tri, Rq(:,1),Rq(:,2),Rq(:,3),LF);
set(ht,'FaceAlpha',faceAlpha)

hold on
qh=quiver3(Rq(:,1),Rq(:,2),Rq(:,3),G(1,:)',G(2,:)',G(3,:)');
hold off
set(qh,'Color','r')

axis equal; axis tight;
shading interp;
colormap(cmap);  colorbar

view(90,90)
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
