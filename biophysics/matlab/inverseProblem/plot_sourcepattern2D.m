function plot_sourcepattern2D(sourceLoc, q,S)
if nargin < 3,
    S=1;
end

nSources = size(sourceLoc,1);
hh = quiver(sourceLoc(:,1), sourceLoc(:,2), q(1:2:end), q(2:2:end), S); hold on
hhm = plot(sourceLoc(:,1), sourceLoc(:,2),'k.'); hold off

% line_handles = get(hh,'Children');
% arrowhead_line = line_handles(2);
% set(arrowhead_line,'LineWidth',4)
set(hh,'Color','k');
set(hhm,'MarkerSize',10)

xlabel('x-axis')
ylabel('y-axis')