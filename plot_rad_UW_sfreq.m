% radar plots for set frequencies over 1 day
clear all; close all

%set frequency
sfreq=0.1; %for T=10s

%load day
day=40;
load ('/data/geophys/scratch/jn6g09/UWDATA/Event_2013_040/LHZ_processedbm_95/LHZ.2013.40','-mat');
%beam has dimensions of theta,frq,slowness,time

%sort out frequnecy vector
frqs=frq(I);

%find frequency index
[jk freqi]=min(abs(frqs-sfreq));
freq=frqs(freqi);
per=1./freq;

%set up day vector
time=([1:Ntime]./Ntime).*24;

%create gridded cartesian coordinates from polar SL and theta
[XX,RAD]=meshgrid(SL,theta*pi/180);
[X1,Y1]=pol2cart(RAD,XX);

%create plot for each timestep
for iday=Ntime:-1:1%[3:3:Ntime-2]
tre=double(squeeze(mean(beam(:,freqi+[-4:1:4],:,iday),2))); %iday+[-2:1:2] %iday
%tre=10*log10(squeeze(mean(tre,3)));
tre=10*log10(squeeze((tre)));
[i,j]=find(tre==max(max(tre)),1);  

figure(iday)%iday-2
%set background colour
set(gcf,'Color',[1,1,1]);
%plot  
h=polar([0 2*pi], [0 max(SL)]); axis('ij'), delete(h), hold on, view([-90 90]); 
pcolor(X1,Y1,real(tre));shading flat;
colorbar; ylabel(colorbar,'power/ dB');
caxis([9 12])
%title(['T= ',num2str(round(per)),'s, day ', num2str(day), '. Time average between ',num2str((round(time(iday-2).*10))./10), ' and ', num2str((round(time(iday+2).*10))./10), 'hr']); 
title(['T= ',num2str(round(per)),'s, day ', num2str(day), '. Time ',num2str((round(time(iday).*10))./10), 'hr']); 
hold on;
[i2, j2]=pol2cart((theta(i)*pi/180),SL(j));
plot(i2,j2,'rs','MarkerSize',15)
text(0.45,-0.6,num2str(tre(i,j)),'FontSize',16,'Color','r')
text(0.35,-0.6,['az= ',num2str(theta(i))],'FontSize',16,'Color','r')
text(0.25,-0.6,['SL= ',num2str(SL(j))],'FontSize',16,'Color','r')
end


        