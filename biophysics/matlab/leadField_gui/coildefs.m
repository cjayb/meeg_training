function [rC, w, er] = coildefs(coilname)

% hack the centerpoint in with weight zero!

if strcmp(coilname, 'mag')
    rC = [-5.25, -5.25, 0.3;...
        -5.25, 5.25, 0.3;...
        5.25, 5.25, 0.3;...
        5.25, -5.25, 0.3] / 1000;
%         5.25, -5.25, 0.3;...
%         0,        0, 0.3] / 1000;
    w = [ones(1,4) / 4];
elseif strcmp(coilname(2:end), 'grad')
    rC = [-8.4, -8.4, 0.3;...
          -8.4,  8.4, 0.3 ;...
           8.4,  8.4, 0.3;...
           8.4, -8.4, 0.3] / 1000;
%            8.4, -8.4, 0.3;...
%              0,    0, 0.3] / 1000;
    if coilname(1) == 'x'
        w = [1, 1, -1, -1] ./ (16.8/1000) ;
    elseif coilname(1) == 'y'
        w = [-1, 1, 1, -1] ./ (16.8/1000) ;
    end
end

er = repmat([0, 0, 1], size(rC,1), 1);
