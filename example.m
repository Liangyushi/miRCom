
%addpath(genpath('/home/amax/amax/lptest/test/'))
load('data.mat');

AUC=cell(35,1);
AUPR=cell(35,1);
testfinal=cell(35,1);
ROC=cell(35,1);
PR=cell(35,1);

%% 5-fold
Pos=find(preX == 1);
pp=size(Pos,1);
index=1:pp;
Indices=crossvalind('Kfold',index,5);
Pos=[Indices Pos];
Pos=[Pos ones(pp,1) zeros(pp,1)];

Neg=find(preX==0);
qq=size(Neg,1);
index=1:qq;
Neg=[index' Neg];
Neg=[Neg zeros(qq,1) zeros(qq,1)];
randnum=randperm(qq);
Nega=Neg(randnum(1:(pp*1)),:);%equal to positive sample

ppa=size(Nega,1);
index=1:ppa;
Indices=crossvalind('Kfold',index,5);
Nega(:,1)=Indices;
clear index Indices Neg

k=1;
Negb=Nega(find(Nega(:,1)==k),:); 
trainX=preX;
A=Pos(find(Pos(:,1)==k),:);
trainX(A(:,2:4))=0;
    
%% Update dm association
dmtest=dm;
dmtest(A(:,[4 2]))=0;
dmtest(A(:,[4 3]))=0;

%% Re Compute Similarity
diseaseSim = constructDiseasesimCell(meshSim,dtSim,dsSim,dmtest, 1, 1);
miRNASim = constructMiRNAsimCell(seqSim,mtSim,meshSim,msSim,dmtest, 1, 1);  

%% Update tensor
Xsparse=sptensor(trainX);
W = trainX;
T=double(trainX);
Omega = (T>0) ;

%% Implement

maxIter= 100;
alpha = 1;        
beta = 0.1;   
K =30;        
R =120;         

%% Obtain Auxiliary features
[LA,Ui] = MyLapalacian(miRNASim,R);
[LB,Vj] = MyLapalacian(diseaseSim,R);

[Xnew,Comp,~,~,~,~] = miRCom(Xsparse,W,miRNASim,diseaseSim,LA,LB,Ui,Vj,maxIter,K,R,alpha,beta);

F=Xnew;
final=[A;Negb];
final(:,6)= F(final(:,2:4));
test=final(:,5:6);
final(:,7)=(mapminmax(final(:,6)', 0, 1))';
testfinal{n}=final;

test=sortrows(test,2);
test(:,3)=(mapminmax(test(:,2)', 0, 1))';
[auc,fpr_tpr,aupr,r_p,~] = Metrics(test(:,1),test(:,3));
fprintf('k:%d : AUC:%.5f ; AUPR:%.5f\n',k,auc,aupr);

AUC{n}=[k auc];
AUPR{n}=[k aupr];
ROC{n}=fpr_tpr;
PR{n}=r_p;
fprintf('Done!\n');

%save('Exampleresults.mat');

