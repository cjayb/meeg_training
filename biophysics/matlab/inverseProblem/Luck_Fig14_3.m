%% draw the scalp surface
figure(101); clf
headR = 0.092;
headTh = linspace(-pi/1.6, pi/1.6, 50);
plot(headR*sin(headTh), headR*cos(headTh), 'Linewidth', 2)
axis equal
hold on
ar = [];
%% generate source space
% kids = get(gca, 'Children');
% delete(kids(1));
ell = struct('majax', [-0.06, 0.01; -0.01, 0.07], ...
             'excent', 0.7, 'theta', [-0.2*pi, 1.2*pi, 16]);
ell(2).majax = [-0.075, -0.01; -0.01, -0.01];
ell(2).excent = 0.8;
ell(2).theta = [pi/2, 2*pi, 14];

sourcePos = [];
sourceOri = [];
for ii=1:length(ell)
    x1 = ell(ii).majax(1,1); x2 = ell(ii).majax(2,1);
    y1 = ell(ii).majax(1,2); y2 = ell(ii).majax(2,2);
    e = ell(ii).excent;
    th1 = ell(ii).theta(1); th2 = ell(ii).theta(2);
    Nth = ell(ii).theta(3);
%     x1 = ell1(1,1); x2 = ell1(2,1);
%     y1 = ell1(1,2); y2 = ell1(2,2);
    
    a = 1/2*sqrt((x2-x1)^2+(y2-y1)^2);
    b = a*sqrt(1-e^2);
    t = linspace(th1,th2, Nth);
    X = a*cos(t);
    Y = b*sin(t);
    w = atan2(y2-y1,x2-x1);
    x = (x1+x2)/2 + X*cos(w) - Y*sin(w);
    y = (y1+y2)/2 + X*sin(w) + Y*cos(w);
    z = zeros(size(x));
    u = cos(t+w); v = sin(t+w); w = zeros(size(u));
    
    sourcePos = [sourcePos; [x(:), y(:), z(:)]];
    sourceOri = [sourceOri; [u(:), v(:), w(:)]];
end
% tweaking: to the the SS like in the figure, swap a couple
% of points!
swaps = [16, 17; 15, 18];
for ii = 1:size(swaps, 1)
    tmp = sourcePos(swaps(ii,1),:);
    sourcePos(swaps(ii,1), :) = sourcePos(swaps(ii,2), :);
    sourcePos(swaps(ii,2), :) = tmp;
end
%% plot source space
plot(sourcePos(:,1),sourcePos(:,2),'m-*', 'Linewidth', 2)
quiver(sourcePos(:,1),sourcePos(:,2),...
       sourceOri(:,1),sourceOri(:,2))

labelInd = [ 0, 5, 10, 15, 16, 20, 25, 29];
labelOff = [-1,-2,  1,  1,  1,  1,  2, -1];
offset = 0.0075;
for ll = 1:length(labelInd)
    if abs(labelOff(ll)) == 1
        offs = [sign(labelOff(ll))*offset, 0];
    else
        offs = [0, sign(labelOff(ll))*offset];
    end
    pos = sourcePos(labelInd(ll)+1, 1:2) + offs; 
    text(pos(1), pos(2), int2str(labelInd(ll)), 'FontSize', 12)
end

%% draw electrodes
elR = headR + 0.0025;
elTh = linspace(-pi/1.7, -0.3, 7)';
elPos = [elR*sin(elTh), elR*cos(elTh), zeros(size(elTh))];
plot(elPos(:,1), elPos(:,2), 'bo', 'MarkerSize', 20)
for ee = 1:size(elPos,1)
    text(elPos(ee, 1), elPos(ee,2), int2str(ee-1), 'FontSize', 12)
end

%% calculate lead field(s)
nSensors = length(elTh);
nSources = size(sourcePos, 1);
% G=zeros(nSensors,nSources);
% for ii = 1:nSensors
%     for jj = 1:nSources
%         % /*     Dipole location    (r): cm                                     */
%         % /*     Dipole moment      (m): mA cm                                  */
%         % /*     Electrode location (s): cm                                     */
%         % /*     *** Output voltage    : mV                                     */
%         r = [100*sourcePos(jj,:), 0];
%         m = [sourceOri(jj,:), 0];
%         s = [100*elPos(ii,:), 0];
%         V = 1000*fastpotential_sun(r, m, s); 
%         G(ii,jj) = V;
%     end
% end
G = calc_potentials(sourcePos, sourceOri, elPos, 1);
%% test null space
[U,S,V] = svd(G);

deepsrc = V(:, 30);  % heuristic
% meas = zeros(nSensors,1);
% for ii = 1:nSensors
%     for jj = 1:nSources
%         % /*     Dipole location    (r): cm                                     */
%         % /*     Dipole moment      (m): mA cm                                  */
%         % /*     Electrode location (s): cm                                     */
%         % /*     *** Output voltage    : mV                                     */
% %         r = [100*sourcePos(jj,:), 0];
% %         m = deepsrc(jj)*[sourceOri(jj,:), 0];
% %         s = [100*elPos(ii,:), 0];
%         r = 100*sourcePos(jj,:);
%         m = 1e5*deepsrc(jj)*sourceOri(jj,:);
%         s = 100*elPos(ii,:);
%         V = 1000*fastpotential_sun(r, m, s); 
%         meas(ii) = meas(ii) + V;
%     end
% end

sourceMom = repmat(deepsrc, 1, 3).*sourceOri;
meas = calc_potentials(sourcePos, sourceMom, elPos);

%% easy position
dipPos = [-0.05, 0.05, 0];
dipMom = 10e-9 * [-1/sqrt(2), 1/sqrt(2), 0];

%% tricky pos, -x
% with low noise, consistent error
% with high noise, more random error
dipPos = [-0.05, 0.007, 0];
dipMom = 10e-9 * [-1, 0, 0];

%% deep pos
% even with minimal noise, substantially wrong solution has low error
dipPos = [-0.012, -0.010, 0];
dipMom = 10e-9 * [1, 0, 0];

%% plot dipole
figure(101)
if ~isempty(ar), delete(ar); end
ar = arrow(dipPos(1:2)+ dipMom(1:2)/2, dipPos(1:2) + dipMom(1:2));
set(ar,'FaceColor', 'c')
%%
dipoleV = calc_potentials(dipPos, dipMom, elPos);
SNR = 30;
norm_name = 'L2-norm';
% norm_name = 'L1-norm';
% norm_name = 'Tikhonov-regularisation';


noiseV = norm(dipoleV)/SNR * randn(size(dipoleV));
measV = dipoleV + noiseV;

if strcmp(norm_name, 'L2-norm')
    estimX = pinv(G)*measV;  % L2-norm
elseif strcmp(norm_name, 'L1-norm')
    estimX = G\measV;  % L1-norm
elseif strcmp(norm_name, 'Tikhonov-regularisation')
    lambda = 1e5;
    estimX = G'*inv(G*G' + lambda*eye(size(G,1)))*measV;  % Tikhonov
    norm_name = [norm_name, ' \lambda = ', num2str(lambda)];
end

predV = G*estimX;

RV = norm(measV - predV)/norm(measV)*100;

figure(201); clf
% subplot(2,1,1); plot([0:nSensors-1], dipoleV, 'Linewidth', 2)
subplot(2,1,1); plot([0:nSensors-1], 1e6*measV, '-', 'Linewidth', 2); hold on
subplot(2,1,1); plot([0:nSensors-1], 1e6*dipoleV, '--', 'Linewidth', 2)
subplot(2,1,1); plot([0:nSensors-1], 1e6*predV, 'c--', 'Linewidth', 2)
% legstr = sprintf('RV: %.1g\%', RV);
legend('meas', 'clean', 'estim', 'Location', 'best')
ylabel('Electric potential (\muV)')

subplot(2,1,2); plot([0:nSources-1], 1e9*estimX, 'Linewidth', 2); hold on
ylabel('Source estimate (nA m)')
title(norm_name)

%%
figure(201);
newestimX = estimX;
for ss = 8:30,
    newestimX = newestimX + max(estimX)*(2*rand()-1)*V(:, ss);
end
nullpredV = G*newestimX;
subplot(2,1,1); plot([0:nSensors-1], 1e6*nullpredV, 'y*', 'Linewidth', 2)
subplot(2,1,2); plot([0:nSources-1], 1e9*newestimX, 'Linewidth', 2)
