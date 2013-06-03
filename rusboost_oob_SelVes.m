%***************************
%*  RUSBoost:
%*  Out of bag
%   289 features (distance from ostium)
%   
%*  2013,6,2
%*  References :
%*  Seiffert et al, 
%*  
%***************************
%
%******************************

clear all; clc; close all;
mkdir models
%% Flags
%% Training mode 
% 0-> binary 0=healthy,soft|1=calc,mix
% 1-> binary 0=healthy,soft,mix|1=calc
% 2-> binary 0=healthy|1=calc,mix,narr>50
% 3-> multi 0=healthy|1=calc,mix|2=narr>50
% 4-> binary 0=healthy|1=grade narrowing>50
trainingMode=0;


%% file lecture
if trainingMode==0
    trainFile='tm0_oob_L5_Train289'
elseif trainingMode==1
   trainFile='BC_L5_allvesselsTrain_OnlyC'
   %trainFile='BC_Circle_allvesselsTrain_OnlyC'
elseif trainingMode==2
   %trainFile='tm2_L5_Train288'
   trainFile='tm2_L5_Segments288'
   %trainFile='BC_Circle_allvesselsTrain_OnlyC'
elseif trainingMode==3
   trainFile='tm3_L5_Segments288'
   %trainFile='BC_Circle_allvesselsTrain_OnlyC'
elseif trainingMode==4
   trainFile='tm4_L5_Segments288'
   %trainFile='BC_Circle_allvesselsTrain_OnlyC'     
else
    trainFile='MC_L5_allvesselsTrain'
end

load(strcat('tmp_data/',trainFile,'.mat'));

tabulate(yTrain)
tabulate(yTest)


%% RUSBOOST
cltree = ClassificationTree.template('minleaf',5);
tic
rusTree = fitensemble(trainData,yTrain,'RUSBoost',500,cltree,...
    'LearnRate',0.1,'nprint',100,'type','classification');
toc
tic
rusTree2 = fitensemble([trainData; testData],[yTrain;yTest(:,4)] ,'RUSBoost',500,cltree,...
    'LearnRate',0.1,'nprint',100,'type','classification','kfold',4)%,...
toc



%% Test
figure;
plot(loss(rusTree,testData,yTest(:,4),'mode','cumulative'));
grid on;
hold on;
plot(kfoldLoss(rusTree2,'mode','cumulative'),'r.');
hold off;
xlabel('Number of trees');
ylabel('Classification error');
legend('Test','4-fold Cross-validation','Location','NE');

% check confusion matrix
tic
Yfit = predict(rusTree,testData);
toc
tab = tabulate(yTest(:,4));
cm=confusionmat(yTest(:,4),Yfit)
cm2=bsxfun(@rdivide,cm,tab(:,2))*100

% Measures
TN=cm(1,1)
TP=cm(2,2)
FN=cm(1,2)
FP=cm(2,1)

SEN=TP/(TP+FN)
FPR=FP/(FP+TN);
SPC=TN/(FP+TN)
ACC=(TP+TN)/(TP+FN+FP+TN)
PPV=TP/(TP+FP)
NPV=TN/(FN+TN)


save models/rb500_tm0_oob_SelVes_289 rusTree rusTree2 cm cm2


%% Roc curve
% binary class
[fpr,tpr] = perfcurve(yTrain(istest),Yfit,1);
% Multi class
%[fpr,tpr] = perfcurve(yTrain(istest),Yfit,2,'negClass',[0 1 3]);
figure;
plot(fpr,tpr);
xlabel('False Positive Rate');
ylabel('True Positive Rate');

roc(yTrain(istest),Yfit)

plotconfusion(yTrain(istest),Yfit);

