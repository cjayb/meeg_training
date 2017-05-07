function V = calc_field(src_vec, sig, X, Y)

%sig = 1; %conductivity = 1
const = 1/(4*pi*sig);

V = 0; % initialize the value of V to zero

for ii = 1:size(src_vec,1) % loop over all rows (i.e., coordinate pairs)
    XX = X-src_vec(ii,1);    % calculate the distance of each point X, from the x-position of the current source/sink
    YY = Y-src_vec(ii,2);    % same for Y
    
    V = V + const*src_vec(ii,3)./sqrt(XX.^2 + YY.^2);   % update the value of V
end
V = V*1000; % mV
