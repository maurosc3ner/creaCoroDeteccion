clear all; clc; close all;
mkdir models
%% Flags
%% Training mode 
% 0-> binary 0=healthy,soft|1=calc,mix
% 1-> binary 0=healthy,soft,mix|1=calc
trainingMode=3;
%% Partition mode
% 0-> manual
% 1-> cvpartition
partitionMode=1;

%% file lecture
if trainingMode==0
    trainFile='BC_L5_Train288'
    %trainFile='BC_Circle_allvesselsTrain'
elseif trainingMode==1
   trainFile='BC_L5_allvesselsTrain_OnlyC'
   %trainFile='BC_Circle_allvesselsTrain_OnlyC'
elseif trainingMode==2
   trainFile='tm2_L5_Train288'
   %trainFile='BC_Circle_allvesselsTrain_OnlyC'
elseif trainingMode==3
   trainFile='tm3_L5_Train288'
   %trainFile='BC_Circle_allvesselsTrain_OnlyC'   
else
    trainFile='MC_L5_allvesselsTrain'
end

load(strcat('tmp_data/',trainFile,'.mat'));

tabulate(yTrain)
if partitionMode==0
    %% Barajamiento de datos
    % Extracting calc data points
    idx = (yTrain==1);
    calc_data = trainData(idx,:); 
    row_calc = size(calc_data,1);

    % Extracting soft data points
    idx2 = (yTrain==2);
    soft_data = trainData(idx2,:); 
    row_soft = size(soft_data,1);

    % Extracting healthy data points
    idx3 = (yTrain==0);
    heal_data = trainData(idx3,:);
    row_heal = size(heal_data,1);

    % Random permuation of positive and negative data points
    c = randperm(row_calc);
    s = randperm(row_soft);
    h = randperm(row_heal);

    % 66-33 split for training and test
    tstsf = s(1:round(row_soft/5));
    tstcf = c(1:round(row_calc/5));
    tsthf = h(1:round(row_heal/5));
    trsf = setdiff(s, tstsf);
    trcf = setdiff(c, tstcf);
    trhf = setdiff(h, tsthf);

    %% cvparition for rustree 1
    train_data = [calc_data(trcf,:);...
                    soft_data(trsf,:);...
                    heal_data(trhf,:)];
    y_train_data = [ ones(numel(trcf),1); ...
                        ones(numel(trsf),1)*2;
                        zeros(numel(trhf),1) ];

    test_data = [calc_data(tstcf,:);...
                    soft_data(tstsf,:);...
                    heal_data(tsthf,:)];
    y_test_data = [ ones(numel(tstcf),1); ...
                        ones(numel(tstsf),1)*2;
                        zeros(numel(tsthf),1) ];

    tabulate(y_train_data)
    tabulate(y_test_data)

    %% cvparition for rustree 2
    train_data2 = [calc_data(:,:);...
                    soft_data(:,:);...
                    heal_data(:,:)];
    y_train_data2 = [ ones(size(calc_data,1),1); ...
                        ones(size(soft_data,1),1)*2;
                        zeros(size(heal_data,1),1) ];

    tabulate(y_train_data2)
    %% RUSBoost 
    cltree = ClassificationTree.template('minleaf',5);
    tic
    rusTree1 = fitensemble(train_data,y_train_data,'RUSBoost',1000,cltree,...
        'LearnRate',0.1,'nprint',100,'type','classification')%,...
    toc
    %% RUSBoost with kfolds=4
    tic
    rusTree2 = fitensemble(train_data2,y_train_data2,'RUSBoost',500,cltree,...
        'LearnRate',0.1,'nprint',100,'type','classification','kfold',4)%,...
    toc
    
elseif partitionMode==1
    part = cvpartition(yTrain,'holdout',0.4);
    istrain = training(part); % data for fitting
    istest = test(part); % data for quality assessment
    tabulate(yTrain(istrain))
    %% RUSBoost 
    cltree = ClassificationTree.template('minleaf',1);
    tic
    rusTree1 = fitensemble(trainData(istrain,:),yTrain(istrain),'RUSBoost',250,cltree,...
        'LearnRate',0.1,'nprint',100,'type','classification')%,...
    toc
    %% RUSBoost with kfolds=4
    tic
    rusTree2 = fitensemble(trainData,yTrain,'RUSBoost',250,cltree,...
        'LearnRate',0.1,'nprint',100,'type','classification','kfold',4)%,...
    toc
end






%% Performance curves
figure;
plot(loss(rusTree1,trainData(istest,:),yTrain(istest),'mode','cumulative'));
%plot(loss(rusTree1,test_data,y_test_data,'mode','cumulative'));
grid on;
hold on;
plot(kfoldLoss(rusTree2,'mode','cumulative'),'r.');
hold off;
xlabel('Number of trees');
ylabel('Classification error');
legend('Test','Cross-validation','Location','NE');

Yfit = predict(rusTree1,trainData(istest,:));
[X,Y]=perfcurve(yTrain(istest),Yfit,[1]);
figure;
plot(X,Y);

[fpr,accu,thre]=perfcurve(yTrain(istest),Yfit,[1],'ycrit','accu');
figure;
plot(thre,accu);

% check confusion matrix
tic
Yfit = predict(rusTree,test_data);
toc
tab = tabulate(y_test_data);
cm=confusionmat(y_test_data,Yfit)
cm2=bsxfun(@rdivide,cm,tab(:,2))*100

save models/rb500_TM3_MC_MP66 rusTree cm cm2

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

