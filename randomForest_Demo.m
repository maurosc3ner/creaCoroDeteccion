%***************************
%*  Random Forest:
%*  Out of bag
%   289 features (distance from ostium)
%   
%*  2013,6,2
%*  References :
%*  , 
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
tic
rf = TreeBagger(250,trainData,yTrain,...
    'method','classification','NVarToSample','all','nprint',100,...
    'minleaf',5,'oobvarimp','on')
toc

%% Test
figure1=figure('Color',[1 1 1]);
plot(oobError(rf));
grid on;
xlabel('Number of grown trees');
ylabel('Classification error');
legend('Random Forest','Location','NE');

% check confusion matrix
tic
Yfit = predict(rf,testData);
toc
tab = tabulate(yTest(:,4));
cm=confusionmat(yTest(:,4),str2num(cell2mat(Yfit)))
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

save models/randomForest rf 

