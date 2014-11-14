clear all
global EV
global U
global epsilon0


ev=dir(['/data/geophys/scratch/jn6g09/CASCADIADATA/Event_2013_018*']); %list of Event directories
ddrr='/data/geophys/scratch/jn6g09/CASCADIADATA/'; %directory containing Event directories

%CHANGE lines 81 and 99 sta.name(1:5) 

for kk=1:length(ev); %loop over Events
    
    clearvars -except EV U epsilon0 ev ddrr kk

stnfile=[ddrr,ev(kk).name,'/LDH_processed1/stations.LDH']
if exist(stnfile) %carry on if LDH is present on this date
load(stnfile, '-mat')



%load -mat /Volumes/JFNexternaldrive/Seismic/CASCADIAresponse/Jan/Testday/Event_2012_266/processed1/stations.LHZ  %change this to pick up the date you want
%figure(1);clf;hold on; plot([infom(:).slon], [infom(:).slat],'r*');grid on
% FOR CAR LOcation is not in file
%if (strmatch(infom(154).staname(1:3),'CAR'))
%infom(154).slat=35.30819;
%infom(154).slon=	-119.84583;
%end
%%

% full grid
%Japan
%LonLref=-180; LonUref= 180;  LatLref= 0; LatUref= 90;  
LonLref=-134; LonUref=-125; LatLref= 0; LatUref= 90;




%LonLref=-125; LonUref= -112;  LatLref= 23; LatUref= 46;  
%
plot([LonLref LonLref LonUref LonUref LonLref],[LatLref LatUref LatUref ...
      LatLref LatLref]);   

stacoord=[ [infom.slon]; [infom.slat]];
%%Find the stations which belong to this grid 
ISTA=find( (stacoord(1,:)>=LonLref ) &  (stacoord(1,:)<=LonUref )& ...
   (stacoord(2,:)>=LatLref ) &  (stacoord(2,:)<=LatUref )  );
size(ISTA)
plot([infom(ISTA).slon], [infom(ISTA).slat],'g*')

meanlat=mean([infom(ISTA).slat]); meanlon=mean([infom(ISTA).slon]);

% readstaresp
 
ic=0;
for ista=ISTA
  ic=ic+1;
 [ran aa b]=dist_wh([meanlat [infom(ista).slat]],[meanlon [infom(ista).slon]]);
 if (aa<0),    aa=360+aa;   end  
 xsta(ic)=ran*sin(aa*pi/180); ysta(ic)=ran*cos(aa*pi/180);
end
coord=[xsta; ysta];
%plot(xsta,ysta,'*');figure; plot([infom(ISTA).slon], [infom(ISTA).slat],'g*')
%%


iyr=str2num(ev(kk).name(7:10))
for imonth=str2num(ev(kk).name(12:14)) %Change julian days here!!!!!!!!!!!!!!  0+[1:4] %:-1:6
%%359 
%% Change the julian date above
  iday1=1;
%Japan
%kcomp={'LHZ' 'BHN' 'BHE'};
%kcomp={'LHZ' 'BHN' 'BHE'};
inpath =   [ddrr,ev(kk).name,'/LDH_processed1/'] ; %change the directory to the sac2matfreq beam
%output directory
    %  sta1=stations{ista};
    for ista=ISTA
      sta1=infom(ista).staname;%1:4
      sta2=cellstr(sta1);
      sta3=char(sta2);
      filename= [inpath,sta3,'.',num2str(iyr),'.', num2str(imonth), '.LDH'];
      if  exist(filename,'file')
	eval(['load -mat ' filename]); %seis1(:,:,ic)=squeeze(fseis(1:end,iday1,:)); 
%	      eval(['load ' filename]);
%	      seis1(:,:,ic)=squeeze(fseis(1:end,:)); 
        break
    end
  end
 Nfreq=size(frq(Imin:Imax),2);
 Nsta=length(ISTA);
 Nframes= max(Nsample(imonth,:));

  seis1=zeros(Nfreq, Nframes,length(ISTA),'single');
    ic=0
    for ista=ISTA
      ic=ic+1;
      %  sta1=stations{ista};
      sta1=infom(ista).staname; %(1:4);
      sta2=cellstr(sta1);
      sta3=char(sta2);
      j1=1;
      %    filename= [inpath sta1 kcomp{j1} num2str(234+iday) '.mat'];
       filename= [inpath sta3,'.',num2str(iyr),'.',  num2str(imonth) '.LDH'];
      if  exist(filename,'file')
	eval(['load  -mat ' filename]); seis1(:,1:size(fseis,3),ic)=squeeze(fseis(1:end,iday1,:)); 
%	      eval(['load ' filename]); seis1(:,:,ic)=squeeze(fseis(1:end,:)); 
 size(fseis,3)
      else
  	display([filename ' did not exist' ])
      end
    end
 %   correct for instrument response
 %    I=Imin:Imax;
%    for ista=ISTA
%      if sum(ista==ISTA1)
%         icc=1;
%       elseif sum(ista==ISTA2)%
%	 icc=2;
%       elseif sum(ista==ISTA3)%
%	 icc=3;
%       else%
%	 icc=0  %then in blow up!
%       end
%       for ii=1:Nfreq%
%	 seis1(ii,:,ista)=seis1(ii,:,ista)/stationresp(I(ii),icc);
%       end
%     end
%%

if size(stacoord,2)>1%don't try and do beamforming if there is only one station of data
    I=Imin:Imax;
    theta = (0:2:360)';
    Nfrq= floor(length(I));
    Nsta=length(coord);
    Ntheta=length(theta);
    %SL=0.0:0.01:0.4;%0.0:0.01:0.4;  %OR 10?
   % SL=0.0:0.001:0.1
  %  SL=0.15:0.01:0.4;
  SL=1:0.5:25;
    thetarad=theta*pi/180;
    projection=[sin(thetarad) cos(thetarad)]*coord;
    timestep=5; %CHANGED FROM 5
    Ntime=floor(Nframes/timestep);
    beam=zeros(Ntheta,Nfrq,length(SL),Ntime,'single');
 %    beamab=zeros(Ntheta,Nfrq,length(SL),Ntime,'single');
  % beam=zeros(Ntheta,Nfrq,length(SL));
    icc=0;
 %keyboard
 for sl=SL   %1000:50:4000
     icc=icc+1
     disp(cputime)
     for ifreq=1:Nfrq
         % this should be more correct 30 nov 8 am      freq=frq(I(ifreq)-1);
         freq=frq(I(ifreq));
         omega=freq*2*pi;
         rep=exp((j*omega*sl/1000)*projection)/sqrt(Nsta);
         for itime=1:Ntime
             itc=(itime-1)*timestep+1;
%change in here to do rotation 

             vhelp= squeeze(seis1(ifreq,itc:itc+timestep,:));
             beam(:,ifreq,icc,itime)= sum(abs(rep*vhelp').^2,2);
       %      C=vhelp'*vhelp;
             %           beam(:,ifreq,icc,itime)= sum(abs(rep*vhelp').^2,2);

             clear EV
%             beamab(:,ifreq,icc,itime)= wncfast(C,rep(:,:)',2,0.2).';
         end
   
     end  %ifreq-v7.3
 end     %SL
%%   
    eval(['save  -v7.3 /data/geophys/scratch/jn6g09/CASCADIADATA/',ev(kk).name,'/LDH_processedbmIGW/','LDH.',num2str(iyr),'.', num2str(imonth) ' beam  frq I  imonth timestep Ntime Nfrq SL theta Ntheta'])
%%
end%end if number of stations>1 loop
end%end imonth loops (loop not needed?)
end %end if file exists loop
end %end Event loop

