clear all
close all

spectr=0;
spectrplot=0;
seisplot=1;
st='G17B';
yr=2013;
type='LDH'

%cd /Volumes/JFNexternaldrive/Seismic/scripts_use/Spectra/
%directory containing event directories
dr='/data/geophys/scratch/jn6g09/CASCADIADATA/'


%get Event list
icnt=0
for i=1:31%329:344%32:41
    icnt=icnt+1;
ev(icnt)=dir([dr,'Event_' char(num2str(yr)),'_',char(num2str(i,'%03.f')),'*']); %list of Event directories
end

%Read in data for station st
for i=1:length(ev)
    %read in data for station
if exist([dr,ev(i).name,'/',num2str(type),'_rsp/',st,'.',char(num2str(yr)),'.',ev(i).name(1,12:14),'.',num2str(type),'.SAC.rsp'])
a(:,:,i)=rsac([dr,ev(i).name,'/',num2str(type),'_rsp/',st,'.',char(num2str(yr)),'.',ev(i).name(1,12:14),'.',num2str(type),'.SAC.rsp']);
%create common time axis
[npts,B,NZJDAY,NZYEAR,NZHOUR,NZMIN,NZSEC,NZMSEC, dt] = lh(a(:,:,i),'NPTS','B','NZJDAY','NZYEAR','NZHOUR','NZMIN','NZSEC','NZMSEC','DELTA');
a(1:86400,1,i)=[0:1:86399]+(NZJDAY-1)*24*3600;
else
if i>1
    a(1:86400,1,i)=a(1:86400,1,i-1)+86400;
    a(1:86400,2,i)=nan;
    a(1:86400,3,i)=nan;
elseif i==1
    a(1:86400,1,i)=nan;
    a(1:86400,2,i)=nan;
    a(1:86400,3,i)=nan;
end
end   
end

%concatenate all the days together
span(1:86400,1:3)=a(:,:,1);
for i=2:length(ev)
span=vertcat(span,a(:,:,i));
end
clear date
date=span(:,1)/86400+datenum(yr,1,1,0,0,0);

%plot the seismogram
if seisplot==1
figure(1)
plot(date,span(:,2))
datetick('x','dd/mm','keepticks')
ylabel('velocity (m/s)','FontSize',14)
title(st,'FontSize',14)
set(gca,'FontSize',14)
end

%filter % between 10 and 20 seconds
Ts=1;
Fe=1/Ts;
freq_int=[0.1 0.2];  
[BB,AA]=butter(4,[freq_int]/Fe*2);
span2= filtfilt(BB,AA,detrend(span(:,2)));


%plot the seismogram
if seisplot==1
figure(2)
plot(date,span2)
datetick('x','dd/mm','keepticks')
ylabel('velocity (m/s)','FontSize',14)
title(st,'FontSize',14)
set(gca,'FontSize',14)
end


% stdev=nanstd(span(:,2));
% span2= span(:,:);
% Threshold_STD=1;
% Threshold_Temp= Threshold_STD*stdev;
% II=find(abs(span2(:,2))>Threshold_Temp);
% span2(II,2)=Threshold_Temp*sign(span2(II,2));
% meanfac=nanmean(span2(:,2)); 
% span2(:,2)=span2(:,2)-ones(size(span2(:,2)))*meanfac;
% span(:,2)=span2(:,2);

dv=datevec(date);

if spectr==1

[S F T P]= spectrogram(span(:,2),2^12,[],2^13,1,'yaxis');
T2=T/86400+date(1);

%plot spectrogram
if spectrplot==1
figure(3)
surf(T2,F,10*log10(abs(P)));
%pcolor(T2,F,abs(P));shading flat;
axis tight;
view(0,90);
shading flat
datetick('x','dd/mm','keeplimits','keepticks')
title(st,'FontSize',14)
ylabel('frequency (Hz)','FontSize',14)
colorbar
%caxis([-40 60])
set(gca,'FontSize',14)
end




seis_date=date;
seis_dv=dv;
seis_ts=span;
seis_spec_date=T2;
seis_spec_psd=P;
seis_spec_f=F;
seis_spec_s=S;

clearvars B NZ* date dv T2 P F S npts i dt span a T icnt

%Plot energy at 10 seconds
ind=find(seis_spec_f>0.0475*2 & seis_spec_f<0.0525*2);%about 10 seconds
ts2=10*log10(abs(nanmean(seis_spec_psd(ind,:))));
ts3=nanmean(seis_spec_psd(ind,:));

ind2=find(seis_spec_f>0.12 & seis_spec_f<0.13); %about 8 seconds
ts5=nanmean(seis_spec_psd(ind2,:));

ind3=find(seis_spec_f>0.0625 & seis_spec_f<0.0700); %about 15s
ts7=nanmean(seis_spec_psd(ind3,:));

ind4=find(seis_spec_f>0.049 & seis_spec_f<0.051); %about 15s
ts9=nanmean(seis_spec_psd(ind4,:));

figure(4)
plot(seis_spec_date,ts5); %choose ts2 or ts3
xlabel('date','FontSize',14)
ylabel('Energy ((m/s)^2/Hz)','FontSize',14)
title(['Energy at 8s, station ' st],'FontSize',14)
datetick('x','dd/mm','keeplimits','keepticks')
set(gca,'FontSize',14)
ylim([0 0.15e-9])
end
