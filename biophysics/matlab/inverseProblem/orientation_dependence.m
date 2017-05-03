%% evenly space orientation vectors of length 10
nOrients = 3;
Q_nAm = zeros(nOrients,3);
[Q_nAm(:,1),Q_nAm(:,2),Q_nAm(:,3)] = sph2cart(zeros(nOrients,1), linspace(pi/8,pi/2,nOrients)', 10);
r0_mm = [0,0,60]; %
%%
figure(502); clf
colormap(blue_white_red)
set(gcf,'Position',[136    61   666   645])
clear hh msp esp

for qq=1:size(Q_nAm,1)
    esp(qq) = subplot(nOrients,3,(qq-1)*3 + 1); 
    hh(qq) = quiver(0,60,Q_nAm(qq,1), Q_nAm(qq,3));
    set(gca,'Xlim',[-20,20]);set(gca,'Ylim',[50,75]); axis equal
    xlabel('x-coordinate'); ylabel('z-coordinate');
    
    subplot(nOrients,3,(qq-1)*3 + 2); 
    fasteeg(r0_mm/1000, Q_nAm(qq,:)*1e-9);
    if qq==1
        eeg_clim=get(gca,'Clim');
    end
    title(sprintf('EEG: %.1f to %.1f uV', 1e6*eeg_clim(1), 1e6*eeg_clim(2)));


    msp(qq) = subplot(nOrients,3,(qq-1)*3 + 3); 
    fastmeg(r0_mm/1000, Q_nAm(qq,:)*1e-9);
    if qq==1
        meg_clim=get(gca,'Clim');
    end
    title(sprintf('MEG: %.1f to %.1f pT', 1e12*meg_clim(1), 1e12*meg_clim(2)));
end
set(hh,'AutoScaleFactor',1.0)
set(hh,'LineWidth', 2)
set(hh,'Marker', 'o')

set(msp,'Clim',meg_clim)
set(esp,'Clim',eeg_clim)
