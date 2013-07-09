%***************************
%*  AdaBoostM1:
%*  Out of bag
%   289 features (distance from ostium)
%   
%*  2013,6,2
%*  References :
%*  Freud Y., Schapire R. 
%*  
%***************************
%   Author: Esteban Correa, maurosc3ner@gmail.com
%******************************

clear all; clc; close all;
mkdir models

%% file lecture
trainFile='featureMatrix';
load(strcat('tmp_data/',trainFile,'.mat'));

tabulate(yTrain)
tabulate(yTest(:,4))


%% Randon Forest
cltree = ClassificationTree.template('minleaf',5);
tic
adb = fitensemble(trainData,yTrain,'AdaBoostM1',250,cltree,...
    'LearnRate',0.1,'nprint',100,'type','classification');
toc

load models/rb250_tm0_oob_SelVes_288

%% Test
figure1=figure('Color',[1 1 1]);
hold on;
plot(loss(adb,testData,yTest(:,4),'mode','cumulative'),'r-');
grid on;

% check confusion matrix
tic
Yfit = predict(adb,testData);
toc
tab = tabulate(yTest(:,4));
cm=confusionmat(yTest(:,4),Yfit)
cm2=bsxfun(@rdivide,cm,tab(:,2))*100

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
[SEN SPC ACC FPR PPV NPV]

save models/adaboostm1 adb
