clear all; clc; close all;
mkdir models
trainFile='featureMatrix'
load(strcat('tmp_data/',trainFile,'.mat'));

tabulate(yTrain)
tabulate(yTest(:,4))
%Loading models
load models/modelExamples
%load models/rusboost
%load models/randomForest 

%% AdaBoostM1
% check confusion matrix
tic
[Yfit1 adbscore] = predict(adb,testData);
toc
tab = tabulate(yTest(:,4));
cm=confusionmat(yTest(:,4),Yfit1);
cm2=bsxfun(@rdivide,cm,tab(:,2))*100;

% Measures
TN=cm(2,2);
TP=cm(1,1);
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
[SEN SPC ACC PPV NPV]

%% RUSBoost
% check confusion matrix
tic
[Yfit2 rtscore] = predict(rusTree,testData);
toc
tab = tabulate(yTest(:,4));
cm=confusionmat(yTest(:,4),Yfit2);
cm2=bsxfun(@rdivide,cm,tab(:,2))*100;

% Measures
TN=cm(2,2);
TP=cm(1,1);
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
[SEN SPC ACC PPV NPV]

%% Random Forest
% check confusion matrix
tic
[Yfit3 rfscore] = predict(rf,testData);
toc
tab = tabulate(yTest(:,4));
cm=confusionmat(yTest(:,4),str2num(cell2mat(Yfit3)));
cm2=bsxfun(@rdivide,cm,tab(:,2))*100;

% Measures
TN=cm(2,2);
TP=cm(1,1);
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
[SEN SPC ACC PPV NPV]

%% Performance
figure1=figure('Color',[1 1 1]);
hold on;
plot(loss(adb,testData,yTest(:,4),'mode','cumulative'),'r-');
hold on;
plot(loss(rusTree,testData,yTest(:,4),'mode','cumulative'),'k-.');
hold on;
plot(oobError(rf),'b--');
hold off;
title('Methods Comparison');
xlabel('Number of trees');
ylabel('Classification error');
legend('AdaboostM1','RUSBoost','Random Forest','Location','NE');
grid on;
hold on
axis([0 250 0.02 0.2])

%% ROC curve
figure2=figure('Color',[1 1 1]);
[fpr1,tpr1] = perfcurve(yTest(:,4),adbscore(:,1),'0');
plot(fpr1,tpr1,'r-');
[fpr2,tpr2] = perfcurve(yTest(:,4),rtscore(:,1),'0');
hold on;
plot(fpr2,tpr2,'k-.');

[fpr3,tpr3] = perfcurve(yTest(:,4),rfscore(:,1),'0');
hold on;
plot(fpr3,tpr3,'b--');
title('ROC curve with distance');
xlabel('False Positive Rate');
ylabel('True Positive Rate');
legend('AdaboostM1','RUSBoost','Random Forest','Location','NE');
grid on;