function [cov]=covAnalysis(covTA,par,cov,statePSG,nev,binranges)
unk=cov.unk;
nor=cov.nor;
evn=cov.evn;
NOR=par.st.nor;
OSA=par.st.osa;
CSA=par.st.csa;
MSA=par.st.msa;
HYP=par.st.hyp;
osaRg=(evn==OSA) & (unk==0);
csaRg=(evn==CSA) & (unk==0);
norRg=(nor==1) & (unk==0);
msaRg=(evn==MSA) & (unk==0);
hypRg=(evn==HYP) & (unk==0);
unkRg= (unk==1);
osaTAcov=covTA(osaRg);
csaTAcov=covTA(csaRg);
msaTAcov=covTA(msaRg);
norTAcov=covTA(norRg);
hypTAcov=covTA(hypRg);
unkTAcov=covTA(unkRg);

osaTAcov=osaTAcov(isnan(osaTAcov)==0);
msaTAcov=msaTAcov(isnan(msaTAcov)==0);
norTAcov=norTAcov(isnan(norTAcov)==0);
csaTAcov=csaTAcov(isnan(csaTAcov)==0);
hypTAcov=hypTAcov(isnan(hypTAcov)==0);
unkTAcov=unkTAcov(isnan(unkTAcov)==0);
cov.c.osa=osaTAcov;
cov.c.csa=csaTAcov;
cov.c.nor=norTAcov;
cov.c.hyp=hypTAcov;
cov.c.msa=msaTAcov;

binLen=length(binranges);
[hisOTAglb] = histc(osaTAcov,binranges);mOTA=mean(osaTAcov);sOTA=std(osaTAcov);
[hisCTAglb] = histc(csaTAcov,binranges);mCTA=mean(csaTAcov);sCTA=std(csaTAcov);
[hisMTAglb] = histc(msaTAcov,binranges);mMTA=mean(msaTAcov);sMTA=std(msaTAcov);
[hisNTAglb] = histc(norTAcov,binranges);mNTA=mean(norTAcov);sNTA=std(norTAcov);
[hisHTAglb] = histc(hypTAcov,binranges);mHTA=mean(hypTAcov);sHTA=std(hypTAcov);
[hisUTAglb] = histc(unkTAcov,binranges);%mOTA=mean(osaTAcov);sOTA=std(osaTAcov);

end