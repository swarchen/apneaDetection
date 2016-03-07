function [NSAidx,svmStruct]=doSVM_Ar_NSA_artiGen_v4(rt,rf,par,trianSize)
P_NUM=par.name;
osaTHOrt=rt.t.osa;
csaTHOrt=rt.t.csa;
msaTHOrt=rt.t.msa;
norTHOrt=rt.t.nor;
hypTHOrt=rt.t.hyp;
osaABDrt=rt.a.osa;
csaABDrt=rt.a.csa;
msaABDrt=rt.a.msa;
norABDrt=rt.a.nor;
% hypABDrt=rt.a.hyp;

osaTHOrf=rf.t.osa;
csaTHOrf=rf.t.csa;
msaTHOrf=rf.t.msa;
norTHOrf=rf.t.nor;
% hypTHOrf=rf.t.hyp;
osaABDrf=rf.a.osa;
csaABDrf=rf.a.csa;
msaABDrf=rf.a.msa;
norABDrf=rf.a.nor;

oABDrt=[osaABDrt;msaABDrt;csaABDrt];
oTHOrt=[osaTHOrt;msaTHOrt;csaTHOrt];
oABDrf=[osaABDrf;msaABDrf;csaABDrf];
oTHOrf=[osaTHOrf;msaTHOrf;csaTHOrf];
nABDrt=[norABDrt;];
nTHOrt=[norTHOrt;];
nABDrf=norABDrf;
nTHOrf=norTHOrf;

numO=length(osaTHOrt);
numN=length(norTHOrt);
numC=length(csaTHOrt);
numM=length(msaTHOrt);
numH=length(hypTHOrt);
numS=trianSize;
NRT=[0 0];
ORT=[0 0];
numEV=numO+numC+numM;
% numNH=numN;
if numEV<numN    
    moTHO=mean([osaTHOrt;msaTHOrt]);
    soTHO=std([osaTHOrt;msaTHOrt]);
    moABD=mean([osaABDrt;msaABDrt]);
    soABD=std([osaABDrt;msaABDrt]);
    mcTHO=mean(csaTHOrt);
    scTHO=std(csaTHOrt);
    mcABD=mean(csaABDrt);
    scABD=std(csaABDrt);
    diff=numN-numEV;
    fracO=fix(diff*(numO/(numO+numC)));
    fracC=fix(diff*(numC/(numO+numC)));
    %% AR
    oTHOgen=oTHOrt( fix(rand(fracO,1)*(numO+numM-1)+1),1 )+randn(fix(diff*(numO/(numO+numC) ) ),1 )*soTHO/10;
    oABDgen=oABDrt( fix(rand(fracO,1)*(numO+numM-1)+1),1 )+randn(fix(diff*(numO/(numO+numC) ) ),1 )*soABD/10;
    cTHOgen=csaTHOrt( fix(rand(fracC,1)*(numC-1)+1),1 )+randn( fix(diff*(numC/(numO+numC) ) ),1 )*scTHO/10;
    cABDgen=csaABDrt( fix(rand(fracC,1)*(numC-1)+1),1 )+randn( fix(diff*(numC/(numO+numC) ) ),1 )*scABD/10;
    oABDrt=[osaABDrt;csaABDrt;msaABDrt;oABDgen;cABDgen];
    oTHOrt=[osaTHOrt;csaTHOrt;msaTHOrt;oTHOgen;cTHOgen];
    %%FR
    soTHO=std([osaTHOrf;msaTHOrf]);
    soABD=std([osaABDrf;msaABDrf]);
    scTHO=std(csaTHOrf);
    scABD=std(csaABDrf);
    oTHOgen=oTHOrf( fix(rand(fracO,1)*(numO+numM-1)+1),1 )+randn(fix(diff*(numO/(numO+numC) ) ),1 )*soTHO/10;
    oABDgen=oABDrf( fix(rand(fracO,1)*(numO+numM-1)+1),1 )+randn(fix(diff*(numO/(numO+numC) ) ),1 )*soABD/10;
    cTHOgen=csaTHOrf( fix(rand(fracC,1)*(numC-1)+1),1 )+randn( fix(diff*(numC/(numO+numC) ) ),1 )*scTHO/10;
    cABDgen=csaABDrf( fix(rand(fracC,1)*(numC-1)+1),1 )+randn( fix(diff*(numC/(numO+numC) ) ),1 )*scABD/10;
    oABDrf=[osaABDrf;csaABDrf;msaABDrf;oABDgen;cABDgen];
    oTHOrf=[osaTHOrf;csaTHOrf;msaTHOrf;oTHOgen;cTHOgen];
elseif numN<numEV
    mnTHO=mean(norTHOrt);
    snTHO=std(norTHOrt);
    mnABD=mean(norABDrt);
    snABD=std(norABDrt);
    diff=numEV-numN;
    %% AR
    nTHOgen=norTHOrt( fix(rand(diff,1)*(numN-1)+1),1 )+randn( diff,1 )*snTHO/10;
    nABDgen=norABDrt( fix(rand(diff,1)*(numN-1)+1),1 )+randn( diff,1 )*snABD/10;
    nABDrt=[norABDrt;nABDgen];
    nTHOrt=[norTHOrt;nTHOgen];
    %% FR
    snTHO=std(norTHOrf);
    snABD=std(norABDrf);
    nTHOgen=norTHOrf( fix(rand(diff,1)*(numN-1)+1),1 )+randn( diff,1 )*snTHO/10;
    nABDgen=norABDrf( fix(rand(diff,1)*(numN-1)+1),1 )+randn( diff,1 )*snABD/10;
    nABDrf=[norABDrf;nABDgen];
    nTHOrf=[norTHOrf;nTHOgen];
end
size(nTHOrt)
size(oTHOrt)   
% itera=101;
for ii=1:1
    ort=[oABDrt, oTHOrt, oABDrf, oTHOrf];
    nrt=[nABDrt, nTHOrt, nABDrf, nTHOrf];
    [no,~]=size(ort);
    [nn,~]=size(nrt);

    strArray=java_array('java.lang.String', nn+no);%cell(Onum+Cnum+Nnum,1);
    strArray(1:no) = java.lang.String('sa');
    strArray(no+1:no+nn) = java.lang.String('nh');
    TYPE = cell(strArray);
    group=TYPE;
    xdata=[ort; nrt];
    options.MaxIter=100000;
    tic
    svmStruct = svmtrain(xdata,group,'kernel_function','rbf','Options',options,'autoscale',false);
    toc
end

oABDrt=[osaABDrt;msaABDrt;csaABDrt];
oTHOrt=[osaTHOrt;msaTHOrt;csaTHOrt];
nABDrt=[norABDrt;];
nTHOrt=[norTHOrt;];
oABDrf=[osaABDrf;msaABDrf;csaABDrf];
oTHOrf=[osaTHOrf;msaTHOrf;csaTHOrf];
nABDrf=[norABDrf;];
nTHOrf=[norTHOrf;];
numEV=numO+numM+numC;
OG=[oABDrt, oTHOrt, oABDrf, oTHOrf];
NG=[nABDrt, nTHOrt, nABDrf, nTHOrf];
Sample=[OG;NG];
Group=svmclassify(svmStruct,Sample) ;
isev=strcmp(Group(1:numEV),'sa');
isnor=strcmp(Group(numEV+1:end),'nh');

NSAidx.TP=length(find(isev==1));
NSAidx.FN=length(find(isev==0));
NSAidx.FP=length(find(isnor==0));
NSAidx.TN=length(find(isnor==1));
NSAidx.SenSA=length(find(isev==1))/length(isev);
NSAidx.Acc=(length(find(isev==1))+length(find(isnor==1)))/(length(isev)+length(isnor));

end
