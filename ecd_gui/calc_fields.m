function calc_fields(src_vec, dist, ax)
%
%   src_vec is (x,y,I)
%   dist:   distance in a.u. of observation. The sources and sinks are
%   plotted around +/- 1, so a "closeup" is dist=1, whereas the "far field"
%   is around d=10. 
%   dist can also be a vector of distances, in which case the resulting
%   potential distributions for each distance are plotted separately


xx=linspace(-dist,dist,dist*25);   % dist*25 grid points in x
yy=linspace(-dist,dist,dist*25);   % dist*25 grid points in y
[X,Y]=meshgrid(xx,yy);  % Sets matrix for grid system in x and y

axes(ax)
cla; % clear all children

sig = 1; %conductivity = 1
const = 1/(4*pi*sig);

V = 0; % initialize the value of V to zero

for ii = 1:size(src_vec,1) % loop over all rows (i.e., coordinate pairs)
    XX = X-src_vec(ii,1);    % calculate the distance of each point X, from the x-position of the current source/sink
    YY = Y-src_vec(ii,2);    % same for Y
    
    V = V + const*src_vec(ii,3)./sqrt(XX.^2 + YY.^2);   % update the value of V
end
V = V*1000; % mV

Vmax = max(abs(V(:)));
clim = [-Vmax, Vmax]/(dist^2);
% This sets up "good" contours to visualize with. If you don't understand
% the syntax, take a look at what the variable contains
contsteps = [0.99*min(V(:)), fliplr(clim(1)*logspace(-3,0,10)), clim(2)*logspace(-3,0,10), 0.99*max(V(:))];

contourf(X,Y,V,contsteps, 'HitTest','off');
set(gca,'CLim',clim); % force the color scale to be symmetric
colormap(blue_white_red)
ch = colorbar;
set(ch, 'FontSize',12)
axis image
%title(['Potential of linear dipole, viewed at distance=' num2str(dist)])
