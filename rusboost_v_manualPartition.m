clear all; clc; close all;
mkdir models
%% Flags
%% Training mode 
% 0-> binary 0=healthy,soft|1=calc,mix
% 1-> binary 0=healthy,soft,mix|1=calc
% 2-> binary 0=healthy|1=calc,mix,narr>50
trainingMode=2;

%% file lecture
if trainingMode==0
    trainFile='BC_L5_Train288'
    %trainFile='BC_Circle_allvesselsTrain'
elseif trainingMode==1
   trainFile='BC_L5_allvesselsTrain_OnlyC'
   %trainFile='BC_Circle_allvesselsTrain_OnlyC'
elseif trainingMode==2
   %trainFile='tm2_L5_Train288'
   trainFile='tm2_oob_L5_Train288'
   %trainFile='BC_Circle_allvesselsTrain_OnlyC'
elseif trainingMode==3
   trainFile='tm3_L5_Train288'
   %trainFile='BC_Circle_allvesselsTrain_OnlyC'   
else
    trainFile='MC_L5_allvesselsTrain'
end

load(strcat('tmp_data/',trainFile,'.mat'));

tabulate(yTrain)

%% Barajamiento de datos
% Extracting positive data points
idx = (yTrain==1);
pos_data = trainData(idx,:); 
row_pos = size(pos_data,1);

% Extracting negative data points
neg_data = trainData(~idx,:);
row_neg = size(neg_data,1);
  
% Random permuation of positive and negative data points
p = randperm(row_pos);
n = randperm(row_neg);

% 80-20 split for training and test
tstpf = p(1:round(row_pos/3));
tstnf = n(1:round(row_neg/3));
trpf = setdiff(p, tstpf);
trnf = setdiff(n, tstnf);

%probar barajamiento en cylinderModel 
train_data = [pos_data(trpf,:);neg_data(trnf,:)];
y_train_data = [ ones(numel(trpf),1); ...
                    zeros(numel(trnf),1) ];
test_data = [pos_data(tstpf,:);neg_data(tstnf,:)];
y_test_data = [ ones(numel(tstpf),1); ...
                    zeros(numel(tstnf),1) ];
                
                
tabulate(y_train_data)
y_train_data=y_train_data;
y_test_data=y_test_data;

%% cvparition for rustree 2
train_data2 = [pos_data(:,:);...
    neg_data(:,:)];
y_train_data2 = [ ones(size(pos_data,1),1);...
    zeros(size(neg_data,1),1) ];

% tabulate(y_train_data)
%% RUSBOOST
cltree = ClassificationTree.template('minleaf',1);
tic
rusTree = fitensemble(train_data,y_train_data,'RUSBoost',250,cltree,...
    'LearnRate',0.1,'nprint',100)%,...
toc

tic
rusTree2 = fitensemble(train_data2,y_train_data2,'RUSBoost',250,cltree,...
    'LearnRate',0.1,'nprint',100,'type','classification','kfold',4)%,...
toc


%% Test
figure;
plot(loss(rusTree,test_data,y_test_data,'mode','cumulative'));
grid on;
xlabel('Number of trees');
hold on;
plot(kfoldLoss(rusTree2,'mode','cumulative'),'r.');
hold off;
xlabel('Number of trees');
ylabel('Classification error');
legend('Test','Cross-validation','Location','NE');

% check confusion matrix
tic
Yfit = predict(rusTree,test_data);
toc
tab = tabulate(y_test_data);
cm=confusionmat(y_test_data,Yfit)
cm2=bsxfun(@rdivide,cm,tab(:,2))*100

save models/rb250_TM2_oob_MP66 rusTree cm cm2

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