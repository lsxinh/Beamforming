clear all
%this is to read output from the beamfiles

fn=dir(['TESTOUTP/*.beam']);

%143 is 18th Jan

for i=1:length(fn)
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
beam_all(i,:,:,:)=beam;
end
figure(1)
 pcolor(tax,1./(freq/2/pi),log10(pwr)*10);shading flat;ylabel('Period s');
 datetick('x',12);title('Beam Power dB');colormap('jet');colorbar
figure(2)
 pcolor(tax,1./(freq/2/pi),az);shading flat;ylabel('Period s');
 datetick('x',12);title('Back Azimuth ^o');colormap('jet');colorbar
figure(3)
 pcolor(tax,1./(freq/2/pi),1000./sls);shading flat;ylabel('Period s');
 datetick('x',12);title('Velocity m/s');colormap('jet');colorbar



%plot timeseries of azimuth, power at chosen frequency
figure(4)
set(gcf,'Color','w')
subplot(2,1,2)
az2=az(15,:); %15 for 150 seconds
plot(tax,az2)
datetick('x','mm/yy','keeplimits','keepticks')
title('Azimuth at 150s','FontSize',12)
ylabel('Azimuth (\circ)','FontSize',12)
xlabel('Date','FontSize',12)
set(gca,'FontSize',12)
hold on;

subplot(2,1,1)
pwr2=pwr(15,:); %15 for 150 seconds
plot(tax,pwr2)
datetick('x','mm/yy','keeplimits','keepticks')
title('Power at 150s','FontSize',12)
ylabel('Beam Power (dB)','FontSize',12)
xlabel('Date','FontSize',12)
set(gca,'FontSize',12)
hold on;

%add lines to events
subplot(2,1,2)
%line([tax(3) tax(3)],[0 400],'Color','k')
%line([tax(30) tax(30)],[0 400],'Color','k')
%line([tax(67) tax(67)],[0 400],'Color','k')
%line([tax(110) tax(110)],[0 400],'Color','k')
%line([tax(194) tax(194)],[0 400],'Color','k')
%line([tax(259) tax(259)],[0 400],'Color','k')
plot(tax(3),250,'+k','LineWidth',1);text(tax(3),275,'1','FontSize',12);
plot(tax(30),350,'+k','LineWidth',1);text(tax(30),375,'2','FontSize',12);
plot(tax(67),350,'+k','LineWidth',1);text(tax(67),375,'3','FontSize',12);
plot(tax(110),100,'+k','LineWidth',1);text(tax(110),125,'4','FontSize',12);
plot(tax(194),300,'+k','LineWidth',1);text(tax(194),325,'5','FontSize',12);
plot(tax(259),250,'+k','LineWidth',1);text(tax(259),275,'6','FontSize',12);


subplot(2,1,1)
%line([tax(3) tax(3)],[2e6 4.5e6],'Color','k')
%line([tax(30) tax(30)],[2e6 4.5e6],'Color','k')
%line([tax(67) tax(67)],[2e6 4.5e6],'Color','k')
%line([tax(110) tax(110)],[2e6 4.5e6],'Color','k')
%line([tax(194) tax(194)],[2e6 4.5e6],'Color','k')
%line([tax(260) tax(260)],[2e6 4.5e6],'Color','k')
plot(tax(3),4.2e6,'+k','LineWidth',1);text(tax(3),4.4e6,'1','FontSize',12);
plot(tax(30),3.5e6,'+k','LineWidth',1);text(tax(30),3.7e6,'2','FontSize',12);
plot(tax(67),3.5e6,'+k','LineWidth',1);text(tax(67),3.7e6,'3','FontSize',12);
plot(tax(110),3.5e6,'+k','LineWidth',1);text(tax(110),3.7e6,'4','FontSize',12);
plot(tax(194),4.2e6,'+k','LineWidth',1);text(tax(194),4.4e6,'5','FontSize',12);
plot(tax(259),3.5e6,'+k','LineWidth',1);text(tax(259),3.7e6,'6','FontSize',12);

return

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

%for all periods/frequencies on chosen day, get azimuth

%index 9 for 100s
%index 13 for 125s
%index 14 for 135s
%index 15 for 150s
%index 18 for 175s
%index 20 for 200s
dv=datevec(tax);

t=datenum(2012,11,28);
t=find(tax==t)

az_per(1)=az(9,t);
az_per(2)=az(13,t);
az_per(3)=az(14,t);
az_per(4)=az(15,t);
az_per(5)=az(18,t);
az_per(6)=az(20,t);
az_per=az_per';

%plot/get power at chosen direction and chosen frequency for each day
%load dispersion curve
load('disp.mat')
%use curve to find slowness at chosen frequency
per=200; %%CHANGE%%%%%%%%%
%[jk,mi]=min(abs(1./(freq/2/pi)-per));
[jk,mi]=min(abs(period-per));
v=vel(mi);
sl=1000./v;
%find energy at each azimuth at that slowness value
[jk,mi]=min(abs(SL-sl));

shoreward=beam_all(:,mi,136,20);%for 270 degrees X secs %%CHANGE%%%%%%%%
seaward=beam_all(:,mi,46,20); %for 90 degrees X secs  %%CHANGE%%%%%%%%%%%%%
r2=seaward./shoreward;
r=sqrt(seaward./shoreward);

t=datenum(2012,12,9); %%CHANGE%%%%%%%%%
t=find(tax==t);
r(t)
