function plot_fieldpattern3D(grid,tri,field, sensorPos)

%clf;
hp=trisurf(tri, grid(:,1),grid(:,2),grid(:,3),field);
axis equal; axis tight; hold on
shading interp;%  colorbar
set(hp,'LineStyle','none')
colormap(blue_white_red)
%colormap(gray)

if nargin > 3,
    hold on
    hh = plot3(sensorPos(:,1),sensorPos(:,2),sensorPos(:,3),'k.');
end
set(hh,'MarkerSize',8);

view(2)
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
