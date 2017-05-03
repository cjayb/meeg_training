function cmap = blue_white_red(resol, clim, symmetric)
%magmap     Colourmap for visualisation of magnetic fields.
%   CMAP = MAGMAP(CLIM) returns a resolx3 matrix for use with the
%   Matlab 'colormap'-command. The argument CLIM is optional:
%   CLIM = [CLIM(1) CLIM(2)], giving the minimum and maximum
%   values of the data in the figure. If no CLIM is given, the
%   limits are taken from the current axis: clim = get(gca, 'CLim');
%
%USAGE: cmap = blue_white_red(clim);
%
%Created by:    CJ Bailey, 22/05/02. (magmap.m)
%Modified by:   CJ Bailey, 30/03/13.

if nargin < 1,
    resol = 64;
    clim = get(gca, 'CLim');
    symmetric = 1;
elseif nargin < 2,
    clim = get(gca, 'CLim');
    symmetric = 1;
elseif nargin < 3,
    symmetric = 1;
end

if symmetric,
	ranmax = max(abs(clim));
	clim = [-ranmax, ranmax];
end
zero= round(-clim(1)/(clim(2)-clim(1))*resol);

hsvmap = ones(resol,3);
hsvmap(1:zero, 1) = 2/3;
hsvmap(zero+1:end, 1) = 1;
hsvmap(1:zero,2) = [1:-1/(zero-1):0]';
hsvmap(zero+1:end,2) = [0:1/(resol-zero-1):1]';

%valueramp = max(30, round(resol/2));
valueramp = [zero, resol-zero];
maxdark = 0.4;
hsvmap(1:valueramp(1), 3) = [maxdark:(1-maxdark)/(valueramp(1)-1):1]';
hsvmap((resol-valueramp(2)+1):end, 3) = [1:-(1-maxdark)/(valueramp(2)-1):(maxdark)]';

cmap = hsv2rgb(hsvmap);
if nargin > 2,
    set(gca, 'CLimMode', 'manual');
    set(gca, 'CLim', clim);
end
