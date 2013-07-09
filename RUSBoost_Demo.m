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
%   Author: Esteban Correa, maurosc3ner@gmail.com
%******************************

clear all; clc; close all;
mkdir models

%% file lecture
trainFile='featureMatrix';
load(strcat('tmp_data/',trainFile,'.mat'));

tabulate(yTrain)
tabulate(yTest(:,4))

%% MODEL
cltree = ClassificationTree.template('minleaf',5);
tic
rusTree = fitensemble(trainData,yTrain,'RUSBoost',250,cltree,...
    'LearnRate',0.1,'nprint',100,'type','classification');
toc

%% Test
figure1=figure('Color',[1 1 1]);
hold on;
plot(loss(rusTree,testData,yTest(:,4),'mode','cumulative'),'r.');
grid on;
title('RUSBoost');
xlabel('Number of trees');
ylabel('Classification error');
legend('RUSBoost','Location','NE');

% check confusion matrix
tic
Yfit = predict(rusTree,testData);
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
%Sensit Specifi Accuracy PositiPValue NegPValue
[SEN SPC ACC PPV NPV]

save models/rusboost rt




