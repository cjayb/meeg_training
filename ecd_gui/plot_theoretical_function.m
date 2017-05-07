% This script assumes that the variable 'extrema' exists in the workspace,
% which is automatically created by the 'ecd_gui' when a linear range is
% chosen. The purpose of this script is to compare the "measured" potential
% along a line to a theoretical function of your choosing. The function 
% should be created before running this script and called 'func_obj'. You 
% may want to use a 'lambda function', which is easy to specify, e.g.:
%
% func_obj = @(x) exp(x.^2)
% is a function that returns the exponential of the squared input value x:
% >> func_obj(2)
%    ans =
%       54.5982  % == exp(2^2) == exp(4)
%
% 'extrema' is a struct with fields
% 'x_ex':  the distance from origin of the point with the extreme value
% 'V_ex':  the value of the potential at this point
% 'x_max': how far to extrapolate


x = linspace(extrema.x_ex, extrema.x_max, 50);
xmin = min(x);

hold on
colspec = 'm--';
plot(x, extrema.V_ex/func_obj(xmin) * func_obj(x),colspec)
