function add_circle(R, linestyle, linewidth, col)

hold on 
nTheta = 100;
theta = linspace(0,2*pi,nTheta);
plot(R(1)*cos(theta), R(1)*sin(theta), linestyle, 'Color',col, ...
    'LineWidth', linewidth, 'Hittest', 'off');
hold off