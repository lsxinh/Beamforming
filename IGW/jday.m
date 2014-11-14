%this is to do a table look up of ordinal day of the year and translate to
%dd/mm/yyyy formate

function [dd,mm,yyyy]=jday(odd,yyyy);

mtx=zeros(12,31);
icnt=0;
for i=1:12
    nd=eomday(yyyy,i);
    mtx(i,1:nd)=icnt+[1:nd];
    icnt=icnt+nd;
end
[mm,dd]=find(mtx==odd);
