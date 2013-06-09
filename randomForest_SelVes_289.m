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
    trainFile='tm0_oob_L5_SelVes289'
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
tabulate(yTest(:,4))


%% Randon Forest
tic
rf = TreeBagger(250,trainData,yTrain,...
    'method','classification','NVarToSample','all','nprint',100,...
    'minleaf',5,'oobvarimp','on')
toc

tic
rf2 = TreeBagger(250,[trainData; testData],[yTrain;yTest(:,4)],...
    'method','classification','NVarToSample','all','nprint',100,...
    'minleaf',5,'oobvarimp','on','kfold',4)
toc

load models/rb250_tm0_oob_SelVes_288 rusTree rusTree2 cm cm2


%% Test
figure1=figure('Color',[1 1 1]);
hold on
plot(oobError(rf),'b--');
grid on;
hold on;
plot(kfoldLoss(rf2,'mode','cumulative'),'r.');
hold off;
xlabel('Number of grown trees');
ylabel('Random Forest Classification error');
legend('Test (60:40)','4-fold Cross-validation','Location','NE');




% check confusion matrix
tic
Yfit = predict(rf,testData);
toc
tab = tabulate(yTest(:,4));
cm=confusionmat(yTest(:,4),str2num(cell2mat(Yfit)))
cm2=bsxfun(@rdivide,cm,tab(:,2))*100


% Measures
TN=cm(1,1);
TP=cm(2,2);
FN=cm(1,2);
FP=cm(2,1);

%TrueNegat TruePosi FalseNeg FalsePosit
[TN TP FN FP]

SEN=TP/(TP+FN);
FPR=FP/(FP+TN);
SPC=TN/(FP+TN);
ACC=(TP+TN)/(TP+FN+FP+TN);
PPV=TP/(TP+FP);
NPV=TN/(FN+TN);

%Sensit Specifi Accuracy FalsePRate PositiPValue NegPValue
[SEN SPC ACC FPR PPV NPV]

save models/rbVSrf_250_tm0_oob_SelVes_289 rusTree rf 

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

