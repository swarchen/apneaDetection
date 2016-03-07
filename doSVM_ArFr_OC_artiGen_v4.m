function [Group,OCidx,svmStruct,conf]=doSVM_ArFr_OC_artiGen_v4(rt,rf,par,trianSize,svmStructNSA)
P_NUM=par.name;
osaTHOrt=rt.t.osa;
csaTHOrt=rt.t.csa;
% csaTHOrt=csa.rt.t;
msaTHOrt=rt.t.msa;
norTHOrt=rt.t.nor;
hypTHOrt=rt.t.hyp;
osaABDrt=rt.a.osa;
csaABDrt=rt.a.csa;
% csaABDrt=csa.rt.a;
msaABDrt=rt.a.msa;
norABDrt=rt.a.nor;
hypABDrt=rt.a.hyp;
%%%%%%%%%%%%%%%%%%%%%%%%%
osaTHOrf=rf.t.osa;
csaTHOrf=rf.t.csa;
msaTHOrf=rf.t.msa;
norTHOrf=rf.t.nor;
hypTHOrf=rf.t.hyp;
osaABDrf=rf.a.osa;
csaABDrf=rf.a.csa;
msaABDrf=rf.a.msa;
norABDrf=rf.a.nor;
hypABDrf=rf.a.hyp;

oABDrt=[osaABDrt;msaABDrt;];
oTHOrt=[osaTHOrt;msaTHOrt;];
oABDrf=[osaABDrf;msaABDrf;];
oTHOrf=[osaTHOrf;msaTHOrf;];
cABDrt=[csaABDrt;];
cTHOrt=[csaTHOrt;];
cABDrf=[csaABDrf;];
cTHOrf=[csaTHOrf;];
numO=length(osaTHOrt);
numN=length(norTHOrt);
numC=length(csaTHOrt);
numM=length(msaTHOrt);
numH=length(hypTHOrt);
numS=trianSize;
CRT=[0 0 0 0];
ORT=[0 0 0 0];
numEV=numO+numM;
if numEV<numC
    moTHO=mean([osaTHOrt;msaTHOrt]);
    soTHO=std([osaTHOrt;msaTHOrt]);
    moABD=mean([osaABDrt;msaABDrt]);
    soABD=std([osaABDrt;msaABDrt]);
    diff=numC-numEV;
%     pool=randperm(numEV);
    oTHOgen=oTHOrt( fix(rand(diff,1)*(numO+numM-1)+1),1 )+randn(diff,1 )*soTHO/10;
    oABDgen=oABDrt( fix(rand(diff,1)*(numO+numM-1)+1),1 )+randn(diff,1 )*soABD/10;
    oABDrt=[osaABDrt;msaABDrt;oABDgen];
    oTHOrt=[osaTHOrt;msaTHOrt;oTHOgen];
    %%FR
    moTHO=mean([osaTHOrf;msaTHOrf]);
    soTHO=std([osaTHOrf;msaTHOrf]);
    moABD=mean([osaABDrf;msaABDrf]);
    soABD=std([osaABDrf;msaABDrf]);
    diff=numC-numEV;
    oTHOgen=oTHOrf( fix(rand(diff,1)*(numO+numM-1)+1),1 )+randn(diff,1 )*soTHO/10;
    oABDgen=oABDrf( fix(rand(diff,1)*(numO+numM-1)+1),1 )+randn(diff,1 )*soABD/10;
    oABDrf=[osaABDrf;msaABDrf;oABDgen];
    oTHOrf=[osaTHOrf;msaTHOrf;oTHOgen];
    numEV=length(oTHOrt);
elseif numC<numEV
    mcTHO=mean(csaTHOrt);
    scTHO=std(csaTHOrt);
    mcABD=mean(csaABDrt);
    scABD=std(csaABDrt);
    diff=numEV-numC;
    cTHOgen=cTHOrt( fix(rand(diff,1)*(numC-1)+1),1 )+randn( diff,1 )*scTHO/10;
    cABDgen=cABDrt( fix(rand(diff,1)*(numC-1)+1),1 )+randn( diff,1 )*scABD/10;
    cABDrt=[csaABDrt;cABDgen];
    cTHOrt=[csaTHOrt;cTHOgen];
    %%FR
    mcTHO=mean(csaTHOrf);
    scTHO=std(csaTHOrf);
    mcABD=mean(csaABDrf);
    scABD=std(csaABDrf);
    diff=numEV-numC;
    cTHOgen=cTHOrf( fix(rand(diff,1)*(numC-1)+1),1 )+randn(diff,1 )*scTHO/10;
    cABDgen=cABDrf( fix(rand(diff,1)*(numC-1)+1),1 )+randn(diff,1 )*scABD/10;
    cABDrf=[csaABDrf;cABDgen];
    cTHOrf=[csaTHOrf;cTHOgen];
end
size(oTHOrf)
size(cTHOrf)


for ii=1:1
    ort=[oABDrt oTHOrt oABDrf oTHOrf];
    crt=[cABDrt cTHOrt cABDrf cTHOrf];
    [no,~]=size(ort);
    [nc,~]=size(crt);

    strArray=java_array('java.lang.String', nc+no);%cell(Onum+Cnum+Nnum,1);
    strArray(1:nc) = java.lang.String('csa');
    strArray(nc+1:no+nc) = java.lang.String('osa');
    TYPE = cell(strArray);
    group=TYPE;
    xdata=[crt; ort];
    options.MaxIter=100000;
    tic
    svmStruct = svmtrain(xdata,group,'kernel_function','rbf','Options',options,'autoscale',false);
    toc
end
oABDrt=[osaABDrt;msaABDrt;];
oTHOrt=[osaTHOrt;msaTHOrt;];
oABDrf=[osaABDrf;msaTHOrf;];
oTHOrf=[osaTHOrf;msaABDrf;];
numEV=numO+numM;
OG=[oABDrt oTHOrt oABDrf oTHOrf];
CG=[csaABDrt csaTHOrt csaABDrf csaTHOrf];

Sample=[OG;CG];
Group=svmclassify(svmStruct,Sample) ;
isosa=strcmp(Group(1:numEV),'osa');
iscsa=strcmp(Group(numEV+1:end),'csa');

OCidx.TP=length(find(iscsa==1));
OCidx.FN=length(find(iscsa==0));
OCidx.FP=length(find(isosa==0));
OCidx.TN=length(find(isosa==1));
OCidx.SenO=length(find(isosa==1))/length(isosa);
OCidx.SenC=length(find(iscsa==1))/length(iscsa);
OCidx.Acc=(length(find(isosa==1))+length(find(iscsa==1)))/(length(isosa)+length(iscsa))

oABDrt=[osaABDrt;msaABDrt;];
oTHOrt=[osaTHOrt;msaTHOrt;];
oABDrf=[osaABDrf;msaTHOrf;];
oTHOrf=[osaTHOrf;msaABDrf;];
numEV=numO+numM;
OG=[oABDrt oTHOrt oABDrf oTHOrf];
CG=[csaABDrt csaTHOrt csaABDrf csaTHOrf];
NG=[norABDrt norTHOrt norABDrf norTHOrf];
Sample=[NG;OG;CG];
Group1=svmclassify(svmStructNSA,Sample) ;
Group2=svmclassify(svmStruct,Sample) ;
[len,~]=size(Sample);
isnor=strcmp(Group1(1:numN),'nh');
issa=strcmp(Group1(1:numN),'sa');
isosa=strcmp(Group2(1:numN),'osa');
iscsa=strcmp(Group2(1:numN),'csa');
conf.N2N=length(find(isnor==1));
conf.N2O=length( find((issa==1) & (isosa==1)));
conf.N2C=length( find((issa==1) & (iscsa==1)));
isnor=strcmp(Group1(numN+1:numN+numEV),'nh');
issa=strcmp(Group1(numN+1:numN+numEV),'sa');
isosa=strcmp(Group2(numN+1:numN+numEV),'osa');
iscsa=strcmp(Group2(numN+1:numN+numEV),'csa');
conf.O2N=length(find((isnor==1)));
conf.O2O=length( find( (issa==1) & (isosa==1) ));
conf.O2C=length( find( (issa==1) & (iscsa==1) ));
isnor=strcmp(Group1(numN+numEV+1:end),'nh');
issa=strcmp(Group1(numN+numEV+1:end),'sa');
isosa=strcmp(Group2(numN+numEV+1:end),'osa');
iscsa=strcmp(Group2(numN+numEV+1:end),'csa');
conf.C2N=length(find(isnor==1));
conf.C2O=length( find( (issa==1) & (isosa==1) ) );
conf.C2C=length( find( (issa==1) & (iscsa==1) ) );

end
