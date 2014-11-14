% % %% beamforming output for whole year
% %fl=dir(['/data/geophys/scratch/ep8g10/pbm_MSci/20*']);
% clear
% system('ls /data/geophys/scratch/jn6g09/UWDATA/Event*/LHZ_processedbm_95/* > filelistUW2')
% fl=textread('filelistUW2','%s');
% fl=char(fl);
% 
% for i=1:length(fl) 
%  fl2=cellstr(fl(i,:));
%  fl3=char(fl2);
%  load(fl3,'-mat') %change : or 1:end-1
%  %pick out frequency 
%  freqs=frq(I);
%  fr1=0.12;
%  [c fr1i] = min(abs(freqs-fr1));
%  ic=fr1i;
%  tre_pre=double(squeeze(beam(:,ic,:,:)));
% for iday=1:size(beam,4);
%      
%        %average/smooth over slowness and cut down to one time window
%         tre=squeeze(tre_pre(:,:,iday));
%         tre=squeeze(nanmean(tre(:,5:41),2));
%         tre=10*log10(tre);
%         tre(isnan(tre))=0;
%         
%         beamp(1:length(tre),(i-1)*size(beam,4)+iday)=tre(:,1);
% end
% display(fl(i,:));
% end
% 
% eval(['save /data/geophys/scratch/jn6g09/Beamforming/ts_UW_LHZ_met2_12 ' 'beamp fr1']);
% clearvars -except Ntheta Ntime Nfrq I frq theta SL timestep
% save('var_UW_LHZ_met2_12')

%% contour plot
clear all;close all;
load (['/data/geophys/scratch/jn6g09/Beamforming/ts_UW_LHZ_met2_12.mat']); %output from section 1
load (['/data/geophys/scratch/jn6g09/Beamforming/var_UW_LHZ_met2_12.mat']) % beamforming variables- slowness, frequency


system('ls /data/geophys/scratch/jn6g09/UWDATA/Event*/LHZ_processedbm_95/* > filelistUW2')
fl=textread('filelistUW2','%s');
fl=char(fl);


%choose dates to plot
start=[2013 1 1 0 0 0];
stop=[2013 2 1 0 0 0];


%create a time vector
clear year day time_pre yearstrt
for i=1:size(fl,1)
yeari=str2double(fl(i,end-7:end-4));year(i,:)=yeari;
dayi=str2double(fl(i,end-2:end));day(i,:)=dayi;
end
for j=1:size(fl,1)
yearstrt(j)=datenum(year(j),1,1,0,0,0);
time_pre(j)=yearstrt(j)+day(j)-1;
end 

for i=1:size(fl,1)
    time((i-1)*Ntime+1)=time_pre(i);
    for j=1:Ntime-1
        time((i-1)*Ntime+1+j)=time_pre(i)+(1/Ntime)*j;
    end 
end
time=time';
dv1=datevec(time);


%cut data down into chosen date
startnum=datenum(start);
stopnum=datenum(stop);
starti=find(round(10.*time)==round(10.*startnum),1,'first');
stopi=find(round(10.*time)==round(10.*stopnum),1,'last');
time=time(starti:stopi);
dv=dv1(starti:stopi);
beamp=beamp(:,starti:stopi);


%power time series
figure(1)
pcolor(time,theta,beamp); shading flat;
hold on;
y=zeros(length(time))+min(theta);%plotting dots at each time step somehow
%forces it to plot it correctly
plot(time,y,'o','MarkerSize',1)
ylabel('Azimuth '); 
xlabel('date');
caxis([7 9])
colorbar;
ylabel(colorbar,'power (dB)');
datetick('x','dd/mm/yy','keeplimits','keepticks')
title([num2str(fr1)]);
set(gcf,'Color',[1,1,1]);

