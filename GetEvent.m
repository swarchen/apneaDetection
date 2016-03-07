function [nev,t,du]=GetEvent(P_NUM)
[ndata, text, alldata]=xlsread(['D:\MATLAB\Breath_Project\EventList\' P_NUM '.xls'],P_NUM);
[nrow,ncol]=size(alldata);
nCSA=0;
nOSA=0;
nMSA=0;
nHYP=0;
Dc=nan;
Do=nan;
Dm=nan;
Dh=nan;
[row, col] = find(isnan(ndata(:,3)));
% for ii=1:length(row)
%     if strcmp(text{row(ii)+1},'Central apnea') || strcmp(text{row(ii)+1},'Obstructive apnea') ||strcmp(text{row(ii)+1},'Mixed apnea')
%     tmp=str2num(text{row(ii)+1,6});
%     ndata(row(ii),3)= tmp(1);
%     end
% end
t_osa=0;
t_csa=0;
t_msa=0;
t_hyp=0;
for ii=1:nrow
    if  strcmp(text{ii},'Ãö¿O') 
        Ht(1)=hour(alldata(ii,3));
        Ht(1)=mod(Ht(1)+12,24);
%         if Ht(1)>=21
%             Ht(1)=Ht(1)-12; %%23->11
%         end
        Mt(1)=minute(alldata(ii,3));
        St(1)=second(alldata(ii,3));
        format long
        t_lightoff=datenum([num2str(Ht) ':' num2str(Mt) ':' num2str(St)],'HH:MM:SS');
    end    
end

[ndata, text, alldata]=xlsread(['D:\MATLAB\Breath_Project\mod_EventList2\' P_NUM '.xls'],'Sheet1');
[nrow,ncol]=size(alldata);
nCSA=0;
nOSA=0;
nMSA=0;
nHYP=0;
% pause;
for ii=1:nrow
    if strcmp(text{ii,3},'Central apnea')% && (isnan(ndata(ii-1,3))==0)
        nCSA=nCSA+1;
        Hc(nCSA)=hour(alldata{ii,1})+12;
        if(Hc(nCSA)>=21)
            Hc(nCSA)=Hc(nCSA)-12;
        end
        Mc(nCSA)=minute(alldata{ii,1});
        Sc(nCSA)=second(alldata{ii,1});
%         if isnan(ndata(ii-1,3))==1
%             tmp=str2num(text{ii,6});
%             Dc(nCSA)=tmp(1);
%         else
            Dc(nCSA)=ndata(ii-4,4);
%         end
        t_csa(nCSA)=(datenum([num2str(Hc(nCSA)) ':' num2str(Mc(nCSA)) ':' num2str(Sc(nCSA))],'HH:MM:SS')-t_lightoff)*60*60*24;        
%         t_csa=(t_csa-t_lightoff)*60*60*24;%(sec)
%         if t_csa(nCSA)<0
%              t_csa(nCSA)=3600*12+t_csa(nCSA);
%         end
    end
    if strcmp(text{ii,3},'Obstructive apnea')% && (isnan(ndata(ii-1,3))==0)
        nOSA=nOSA+1;
        Ho(nOSA)=hour(alldata{ii,1})+12;
        if(Ho(nOSA)>=21)
            Ho(nOSA)=Ho(nOSA)-12;
        end
        Mo(nOSA)=minute(alldata{ii,1});
        So(nOSA)=second(alldata{ii,1});
%         if isnan(ndata(ii-1,3))==1
%             tmp=str2num(text{ii,6});
%             Do(nOSA)=tmp(1);
%         else
            Do(nOSA)=ndata(ii-4,4);
%         end
%         Do(nOSA)=ndata(ii-1,3);
% format long
% datenum([num2str(Ho(nOSA)) ':' num2str(Mo(nOSA)) ':' num2str(So(nOSA))],'HH:MM:SS')
        t_osa(nOSA)=(datenum([num2str(Ho(nOSA)) ':' num2str(Mo(nOSA)) ':' num2str(So(nOSA))],'HH:MM:SS')-t_lightoff)*60*60*24; 
%         t_osa=(t_osa-t_lightoff)*60*60*24;%(sec)

%         if t_osa(nOSA)<0
%              t_osa(nOSA)=3600*12+t_osa(nOSA);
%         end
    end
    if strcmp(text{ii,3},'Mixed apnea')% && (isnan(ndata(ii-1,3))==0)
        nMSA=nMSA+1;
        Hm(nMSA)=hour(alldata{ii,1})+12;
        if(Hm(nMSA)>=21)
            Hm(nMSA)=Hm(nMSA)-12;
        end
        Mm(nMSA)=minute(alldata{ii,1});
        Sm(nMSA)=second(alldata{ii,1});
%         if isnan(ndata(ii-1,3))==1
%             tmp=str2num(text{ii,6});
%             Do(nOSA)=tmp(1);
%         else
            Dm(nMSA)=ndata(ii-4,4);
%         end
%         Do(nOSA)=ndata(ii-1,3);
        t_msa(nMSA)=(datenum([num2str(Hm(nMSA)) ':' num2str(Mm(nMSA)) ':' num2str(Sm(nMSA))],'HH:MM:SS')-t_lightoff)*60*60*24; 
%         t_osa=(t_osa-t_lightoff)*60*60*24;%(sec)
%         if t_msa(nMSA)<0
%              t_msa(nMSA)=3600*12+t_msa(nMSA);
%         end
    end    
    if strcmp(text{ii,3},'Hypopnea')% && (isnan(ndata(ii-1,3))==0)
        nHYP=nHYP+1;
        Hh(nHYP)=hour(alldata{ii,1})+12;
        if(Hh(nHYP)>=21)
            Hh(nHYP)=Hh(nHYP)-12;
        end
        Mh(nHYP)=minute(alldata{ii,1});
        Sh(nHYP)=second(alldata{ii,1});
%         if isnan(ndata(ii-1,3))==1
%             tmp=str2num(text{ii,6});
%             Do(nOSA)=tmp(1);
%         else
            Dh(nHYP)=ndata(ii-4,4);
%         end
%         Do(nOSA)=ndata(ii-1,3);
        t_hyp(nHYP)=(datenum([num2str(Hh(nHYP)) ':' num2str(Mh(nHYP)) ':' num2str(Sh(nHYP))],'HH:MM:SS')-t_lightoff)*60*60*24; 
%         t_osa=(t_osa-t_lightoff)*60*60*24;%(sec)
%         if t_hyp(nHYP)<0
%              t_hyp(nHYP)=3600*12+t_hyp(nHYP);
%         end
    end
end
nev.osa=nOSA;
nev.csa=nCSA;
nev.msa=nMSA;
nev.hyp=nHYP;
t.osa=t_osa;
t.csa=t_csa;
t.msa=t_msa;
t.hyp=t_hyp;
du.osa=Do;
du.csa=Dc;
du.msa=Dm;
du.hyp=Dh;

end