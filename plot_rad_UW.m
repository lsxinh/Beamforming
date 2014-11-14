%% radar plots of 1 day
clear all; close all
%% Actual beamforming data
% load day
%load (['/data/geophys/scratch/ep8g10/pbm_MSci/2012_302.mat']);
load ('/data/geophys/scratch/jn6g09/UWDATA/Event_2013_020/LHZ_processedbm_95/LHZ.2013.20','-mat');


% %% data reshaped for radar plot (slowness against azimuth)
% isub=0;
 for ic=[5:5:length(I)-5];  
    for iday=3:Ntime-3;%3:30
        % average/smooth over frequency range 
        tre=double(squeeze(mean(beam(:,ic+[-4:1:4],:,iday+[-2:1:2]),2)));
        %average/smooth over windows
        tre=10*log10(squeeze(mean(tre,3)));
% %        tre=tre-max(max(tre));
%         [i,j]=find(tre==max(max(tre)),1);  
%         isub=isub+1; %islow(isub)=SL(j(1));
% 
% %% Plot
% %make a radar plot of the beamformer output, slowness concentric circles,
% %azimuthal angle
% 
time=((iday)-1)*24/length(Ntime)
% %beamRt=10*log((tre./min(min(tre))));for actual data
% % setup grid
% %vector of slowness in s/km
% SL=0:0.0098:0.4;
% % azimuth in deg
% theta=0:2:360;
% 
 [XX,RAD]=meshgrid(SL,theta*pi/180);
 [X1,Y1]=pol2cart(RAD,XX);
%  
figure(1)
%  
%  %set background colour
 set(gcf,'Color',[1,1,1]);
% 
 h=polar([0 2*pi], [0 max(SL)]); axis('ij'), delete(h), hold on, view([-90 90]); 
 pcolor(X1,Y1,real(tre));shading flat;colorbar;
 per=1./frq(I(ic)) 
 title(['T= ',num2str(per),'s','Time=',num2str(time),' hr']); %for plane wave
 pause(.1)
    end
 end
   
%% for set frequencies and times

return
figure(12)
 %set background colour
set(gcf,'Color',[1,1,1]);

%ic =70 ; % for T = 6s => f=  0.1667Hz
%ic =11 ; % f =0.05Hz, T=20s
%ic = 36;% f = 0.1Hz, T=10s
%ic = 31;
%iday =15 ; % 
%iday=27;
%iday = 17;
iday=6;
time_b = ((iday -2)./Ntime).*24; %start time %changed 33 to Ntime
time_e = ((iday +2)./Ntime).*24; %end time %changed 33 to Ntime

ic=25; %f=0.0508Hz, T=19.7s

[XX,RAD]=meshgrid(SL,theta*pi/180);
 [X1,Y1]=pol2cart(RAD,XX);

    tre=double(squeeze(mean(beam(:,ic+[-4:1:4],:,iday+[-2:1:2]),2)));
    tre=10*log10(squeeze(mean(tre,3)));
    [i,j]=find(tre==max(max(tre)),1);  
  
h=polar([0 2*pi], [0 max(SL)]); axis('ij'), delete(h), hold on, view([-90 90]); 
pcolor(X1,Y1,real(tre));shading flat;
colorbar; ylabel(colorbar,'power/ dB');
per=1./frq(I(ic)) 
title(['T= ',num2str(round(per)),'s, day 183, time average between ',num2str((round(time_b.*10))./10), ' and ', num2str((round(time_e.*10))./10), 'hr']); 
        

        
