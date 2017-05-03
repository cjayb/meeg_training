function plotSphField(func_name, rq, q, R, M, plot_streams, ax)

%clim = [-9.7259    9.7259];

stream_col = [0,0.5,0.5];
Nstreams = 15;

boundary_col = [0.5,0.5,0];

% how many grid points
% N = 256;
N = 128;
xlim = [-1.2, 1.2];
ylim = [-1.2, 1.2];
gpx = round(diff(xlim)/2)*N; % even to avoid origin in set
gpy = round(diff(ylim)/2)*N; % even to avoid origin in set
xx=linspace(xlim(1),xlim(2),gpx);   
yy=linspace(ylim(1),ylim(2),gpy);   

% changing the names here so that x->y and y->z, for easier extension to 3D
[Y,Z]=meshgrid(xx,yy);  % Sets matrix for grid system in x and y

% nRadius = 1;
% nTheta = N;
% % generate a mesh in polar coordinates
% [sphr,theta] = meshgrid(R,linspace(0,2*pi,nTheta));
% % transform from polar to Cartesian coordinates
% Y = [Y, sphr.*cos(theta)];
% Z = [Z, sphr.*sin(theta)];

func_obj = eval(['@', func_name, ';']);
V = func_obj(rq, q, R, M, Y, Z) * 1e9;  % convert to nV!

% tmp = [abs(min(V(:)))^(rq(3)/2), abs(max(V(:)))^(rq(3)/2)];
tmp = [abs(min(V(:)))^(1/2.5), abs(max(V(:)))^(1/2.5)];  % heuristic!
clim = [-max(tmp), max(tmp)];

axes(ax); cla

colormap(blue_white_red);
contsteps = [0.99*min(V(:)), fliplr(clim(1)*logspace(-3,0,10)), clim(2)*logspace(-3,0,10), 0.99*max(V(:))];
contourf(Y,Z,V,contsteps, 'LineStyle','-.', 'HitTest','off'); 
%clim = [contsteps(2), contsteps(end-1)];

if plot_streams,
    [FY, FZ] = gradient(V, 0.05);
    grad   = sqrt(FY.^2 + FZ.^2);
    gradVy = -FY./grad; % Note the minus-sign: with that in place, the arrows
    gradVz = -FZ./grad; % will denote the direction of VOLUME current (J^v = - sigma grad V)

    QAngle = sign(asin(q(2)))*acos(q(3));

    % [sphr,theta] = meshgrid(0.1, linspace(pi/2, 5/2*pi,Nstreams));
    [sphr,theta] = meshgrid(0.1, linspace(0, pi,Nstreams)-QAngle);
    sy = rq(2) + sphr.*cos(theta);
    sz = rq(3)*1.1 + sphr.*sin(theta);
    %[sy,sz] = meshgrid(linspace(0.0,0.4,10), linspace(0.2, 0.6,10));
    sh=streamline(Y,Z,gradVy,gradVz,sy, sz);
    set(sh,'Color',stream_col, 'LineWidth',2)
    
    % h=quiver(Y(1:10:end),Z(1:10:end),...
    %     arrowScale(1:10:end).*gradVy(1:10:end),...
    %     arrowScale(1:10:end).*gradVz(1:10:end),'AutoScale','off'); axis image
    %set(h,'AutoScale','off') % this just scales the arrows a bit...
    % set(h,'Color','r', 'LineWidth',1.)
    
end
circleMask = zeros(size(V));
circleAlpha = ones(size(V));
r = sqrt(Y.^2 + Z.^2);
circleAlpha(r <= R) = 0;


if ~strcmp(func_name, 'infiniteSpace')
      imagesc(xx,yy,circleMask, 'AlphaData', circleAlpha, 'HitTest','off');
end

nTheta = 100;
theta = linspace(0,2*pi,nTheta);
plot(R*cos(theta), R*sin(theta), 'Color',boundary_col, ...
    'LineWidth', 4, 'Hittest', 'off');

hold on; axis xy; axis equal;
set(gca, 'CLim', clim); colorbar
% set(gca, 'XLimMode', 'manual')
% set(gca, 'XLim', xlim);set(gca, 'YLim', ylim);
ar = arrow(rq(2:3)-q(2:3)/10, rq(2:3)+q(2:3)/10);
set(ar, 'FaceColor', 'c', 'EdgeColor','c', 'LineWidth', 2)