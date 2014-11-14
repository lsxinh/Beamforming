clear all
%this is to read output from the beamfiles
addpath(genpath('../../IGWPressureCodes/m_map'))
fn=dir(['TESTOUTP/*.beam']);

%143 is 18th Jan 2013
%103 is 9th Dec 2012
%79 is 15th Nov 2012
%58 is 25th Oct 2012
%30 is 27th Sep 2012
%39 is 6th Oct 2012
%110 is 16th Dec 2012
%6 is 3 sep 2012
%194 is 10 March 2013
%259 is 14 May 2013
%92 is 28 Nov 2012
%113 is 19 Dec 2012
%11th May 2013

for i=219%1:length(fn)
    [dd,mm,yyyy]=jday(str2num(fn(i).name(12:14)),str2num(fn(i).name(7:10)));
    tax(i)=datenum(yyyy,mm,dd);
fid=fopen(['TESTOUTP/',fn(i).name],'r');
a=fread(fid,'single');
fclose(fid)
nf=a(1);nSL=a(2);ntheta=a(3);
asum=3;
freq=a(asum+[1:nf]);
asum=asum+nf;
SL=a(asum+[1:nSL]);
asum=asum+nSL;
theta=a(asum+[1:ntheta]);
asum=asum+ntheta;
beam=zeros(nSL,ntheta,nf);
if i==1
    pwr=zeros(nf,length(fn));
    az=pwr;
    sls=pwr;
end
for k=1:nf
    for ii=1:nSL
        for jj=1:ntheta
            asum=asum+1;
            beam(ii,jj,k)=a(asum);
        end
    end
    pwr(k,i)=max(max(beam(:,:,k)));
    if isnan(pwr(k,i))==1
        az(k,i)=NaN;
        sls(k,i)=NaN;
    else
    [ik,jk]=find(squeeze(beam(:,:,k))==max(max(squeeze(beam(:,:,k)))),1);
    az(k,i)=theta(jk);
    sls(k,i)=SL(ik);
    end
end

end
%figure(1)
 %pcolor(tax,1./(freq/2/pi),log10(pwr)*10);shading flat;ylabel('Period s');
 %datetick('x',12);title('Beam Power dB');colormap('jet');colorbar
%figure(2)
% pcolor(tax,1./(freq/2/pi),az);shading flat;ylabel('Period s');
% datetick('x',12);title('Back Azimuth ^o');colormap('jet');colorbar
%figure(3)
% pcolor(tax,1./(freq/2/pi),1000./sls);shading flat;ylabel('Period s');
% datetick('x',12);title('Velocity m/s');colormap('jet');colorbar
%return


%polar plot for each frequency (need to run top part for chosen day first)
[dd,ss]=meshgrid(theta,SL);
figure(90)
%ii=15 for 150s
for ii=15%1:length(freq)
figure(90)
[X,Y]=pol2cart(dd*pi/180,ss);
h=polar([0 2*pi],[0 (max(ss(:)))]); hold on;axis ij;view([-90 90]);
pcolor(X,Y,(beam(:,:,ii)));shading flat;
title(num2str(2*pi./freq(ii)));%pause


end

%plot dispersion curve

%%
%%create dispersion curve (only need to do once)
%period=1./(freq/2/pi); %period in sec
%vel=1000./nanmean(sls,2); %velocity in m/s
%figure(4)
%plot(period,vel);
%ylabel('Velocity (m/s)')
%xlabel('Period (s)')
%save('disp.mat','period','vel');
%%
%load dispersion curve
load('disp.mat')
%use curve to find slowness at chosen frequency
per=150; %%%%CHANGE
%[jk,mi]=min(abs(1./(freq/2/pi)-per));
[jk,mi]=min(abs(period-per));
v=vel(mi);
sl=1000./v;
%find energy at each azimuth at that slowness value
[jk,mi]=min(abs(SL-sl));
EperAz=beam(mi,:,15); %15 for 150s  %%%CHANGE


%load station data
load('../../IGWPressureCodes/stations.mat');
%calculate array centre
meanlat=mean(stninfo(:,1));
meanlon=mean(stninfo(:,2));
%create grid covering Pacific
lat=[90:-0.2:-90];lon=[100:0.2:290];
[LON,LAT]=meshgrid(lon,lat);
%find azimuth between array centre and grid points
E=referenceEllipsoid('earth');
az=azimuth(meanlat,meanlon,LAT,LON,E);
%project energy per azimuth onto azimuth values
bp=interp1(theta,EperAz,az);

%plot
addpath(genpath('../../IGWPressureCodes/m_map'))
figure(5)
m_proj('Robinson','lat',[-80 80],'lon', [100 290])
%m_coast('patch',[0.4 0.4 0.4],'edgecolor','none');
m_grid('tickdir','out','yaxislocation','left','xaxislocation','bottom'...
    ,'ticklen',.02,'fontsize',12);
hold on;
m_pcolor(LON,LAT,bp);shading flat
m_coast('patch',[0.4 0.4 0.4],'edgecolor','none');
title(datestr(tax(i)),'FontSize',14)
set(gca,'FontSize',14)
set(gcf,'Color','w')
colorbar
caxis([1.5e6 3e6])


%zoomed in plot
figure(6)
m_proj('Mercator','lat',[30 60],'lon', [210 260])
%m_coast('patch',[0.4 0.4 0.4],'edgecolor','none');
m_grid('tickdir','out','yaxislocation','left','xaxislocation','bottom'...
    ,'ticklen',.02,'fontsize',12);
hold on;
m_pcolor(LON,LAT,bp);shading flat
m_coast('patch',[0.4 0.4 0.4],'edgecolor','none');
title(datestr(tax(i)),'FontSize',14)
set(gca,'FontSize',14)
set(gcf,'Color','w')
colorbar
caxis([1.5e6 3e6])

%zoomed in plot
%create grid covering Pacific
lat2=[50:-0.05:30];lon2=[200:0.05:250];
[LON2,LAT2]=meshgrid(lon2,lat2);
%find azimuth between array centre and grid points
E=referenceEllipsoid('earth');
az2=azimuth(meanlat,meanlon,LAT2,LON2,E);
%project energy per azimuth onto azimuth values
bp2=interp1(theta,EperAz,az2);
figure(7)
m_proj('Mercator','lat',[39 46],'lon', [230 240])
%m_coast('patch',[0.4 0.4 0.4],'edgecolor','none');
m_grid('tickdir','out','yaxislocation','left','xaxislocation','bottom'...
    ,'ticklen',.02,'fontsize',12);
hold on;
m_pcolor(LON2,LAT2,bp2);shading flat
m_coast('patch',[0.4 0.4 0.4],'edgecolor','none');
title(datestr(tax(i)),'FontSize',14)
set(gca,'FontSize',14)
set(gcf,'Color','w')
colorbar
caxis([1.5e6 3e6])
