function clim_auto = plot_field(src_vec, dist, sig, ax)
%
%   src_vec is (x,y,I)
%   dist:   distance in a.u. of observation. The sources and sinks are
%   plotted around +/- 1, so a "closeup" is dist=1, whereas the "far field"
%   is around d=10. 
%   dist can also be a vector of distances, in which case the resulting
%   potential distributions for each distance are plotted separately


xx=linspace(-dist,dist,100);   % dist*25 grid points in x
yy=linspace(-dist,dist,100);   % dist*25 grid points in y
[X,Y]=meshgrid(xx,yy);  % Sets matrix for grid system in x and y

axes(ax)
cla; % clear all children

V = calc_field(src_vec, sig, X, Y);

Vmax = max(abs(V(:)));
clim_auto = [-Vmax, Vmax]/(dist^1.5);
% This sets up "good" contours to visualize with. If you don't understand
% the syntax, take a look at what the variable contains
contsteps = [0.95*min(V(:)), fliplr(clim_auto(1)*logspace(-3,0,10)), ...
    clim_auto(2)*logspace(-3,0,10), 0.95*max(V(:))];

contourf(X,Y,V,contsteps, 'HitTest','off');
set(gca,'CLim',clim_auto); % force the color scale to be symmetric
colormap(blue_white_red)
ch = colorbar;
set(ch, 'FontSize',12)
axis image
%title(['Potential of linear dipole, viewed at distance=' num2str(dist)])
hold on

% [foo, argmax] = max(V(:));
% [foo, argmin] = min(V(:));
% plot(X(argmax), Y(argmax), '.', 'MarkerEdgeColor', [0.5, 0.5, 0], 'MarkerSize', 20)
% plot(X(argmin), Y(argmin), '.', 'MarkerEdgeColor', [0.5, 0.5, 0], 'MarkerSize', 20)
