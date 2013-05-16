clear all; clc; close all;
load 'tmp_data/BT_L5_allvesselsTrain'
tabulate(yTrain)

%Particion de datos
part = cvpartition(yTrain,'holdout',0.4);
istrain = training(part); % data for fitting
istest = test(part); % data for quality assessment
tabulate(yTrain(istrain))

%% LibSVM
addpath ../libsvm-3.12/matlab
addpath src

% train_featNorm=(trainData(istrain,:) - repmat(min(trainData(istrain,:),[],1),...
%     size(trainData(istrain,:),1),1))*spdiags(1./(max(trainData(istrain,:),[],1)...
%     -min(trainData(istrain,:),[],1))',0,size(trainData(istrain,:),2),size(trainData(istrain,:),2));
% One line normalization
train_featNorm = normalize_features(trainData(istrain,:), 1);
training_instance_matrix = sparse(train_featNorm);


% test_featNorm = (trainData(istest,:) - repmat(min(trainData(istest,:),[],1),...
%     size(trainData(istest,:),1),1))*spdiags(1./(max(trainData(istest,:),[],1)...
%     -min(trainData(istest,:),[],1))',0,size(trainData(istest,:),2),size(trainData(istest,:),2));
test_featNorm=normalize_features(trainData(istest,:), 1);
test_instance_matrix = sparse(test_featNorm);

%% Linear
%train
model1 = svmtrain(yTrain(istrain,:), training_instance_matrix , '-s 0 -t 0 -g 0.1 -c 1');
%predict
%[unused, unused, scores_test] = svmpredict(testing_label_vector, [(1:size(ker_test,1))' ker_test], model);
[predicted_label1, accuracy1, scores_train1] = svmpredict(yTrain(istrain,:),training_instance_matrix, model1);
[predicted_label1_t, accuracy1_t, scores_test1] = svmpredict(yTrain(istest,:),test_instance_matrix, model1);

tab2=ones(2,1);
tab = tabulate(yTrain(istest,1));
for i=1:size(tab,1)
    tab2(i,1)=tab(i,2);
end
cm=confusionmat(yTrain(istest,1),predicted_label1_t);
cm2=bsxfun(@rdivide,cm,tab2(:))*100

% Measures

TP=cm2(1,1)
TN=cm2(2,2)
FP=cm2(1,2)
FN=cm2(2,1)

TPR=TP/(TP+FN)
FPR=FP/(FP+TN)
ACC=(TP+TN)/(TP+FN+FP+TN)
SPC=1-FPR
PPV=TP/(TP+FP)
NPV=TN/(FN+TN)

%save tmp_data/linear_scores predicted_label1 predicted_label1_t accuracy1 accuracy1_t scores_train1 scores_test1

%% polynomial
model2 = svmtrain(yTrain(istrain,:),training_instance_matrix, '-s 0 -t 1 -d 1 -g 0.001 -c 1');
%train
[predicted_label2, accuracy2, scores_train2] = svmpredict(yTrain(istrain,:),training_instance_matrix, model2);
[predicted_label2_t, accuracy2_t, scores_test2] = svmpredict(yTrain(istest,:),test_instance_matrix, model2);
%predict
%save tmp_data/poly_scores predicted_label2 predicted_label2_t accuracy2 accuracy2_t scores_train2 scores_test2
tab2=ones(2,1);
tab = tabulate(yTrain(istest,1));
for i=1:size(tab,1)
    tab2(i,1)=tab(i,2);
end
cm=confusionmat(yTrain(istest,1),predicted_label2_t);
cm2=bsxfun(@rdivide,cm,tab2(:))*100

% Measures
TP=cm2(1,1)
TN=cm2(2,2)
FP=cm2(1,2)
FN=cm2(2,1)

TPR=TP/(TP+FN)
FPR=FP/(FP+TN)
ACC=(TP+TN)/(TP+FN+FP+TN)
SPC=1-FPR
PPV=TP/(TP+FP)
NPV=TN/(FN+TN)

%% rbf
%train
bestcv=0;
bestc=[];
bestg=[];
for log2c = -1:3,
  for log2g = -4:1,
    cmd = ['-v 5 -s 0 -t 2 -c ', num2str(2^log2c), ' -g ', num2str(2^log2g)];
    %cv = svmtrain(heart_scale_label, heart_scale_inst, cmd);
    cv = svmtrain(yTrain(istrain,:), training_instance_matrix , cmd);
    if (cv >= bestcv),
      bestcv = cv; bestc = 2^log2c; bestg = 2^log2g;
    end
    fprintf('%g %g %g (best c=%g, g=%g, rate=%g)\n', log2c, log2g, cv, bestc, bestg, bestcv);
  end
end

2^bestc
2^bestg
model5 = svmtrain(yTrain(istrain,:), training_instance_matrix , '-s 0 -t 2 -g 1 -c 1');
%predict
[predicted_label5, accuracy5, scores_train5] = svmpredict(yTrain(istrain,:),training_instance_matrix, bestcv);
[predicted_label5_t, accuracy5_t, scores_test5] = svmpredict(yTrain(istest,:),test_instance_matrix, bestcv);
%save tmp_data/rbf_scores predicted_label5 predicted_label5_t accuracy5 accuracy5_t scores_train5 scores_test5
tab2=ones(2,1);
tab = tabulate(yTrain(istest,1));
for i=1:size(tab,1)
    tab2(i,1)=tab(i,2);
end
cm=confusionmat(yTrain(istest,1),predicted_label5_t);
cm2=bsxfun(@rdivide,cm,tab2(:))*100

% Measures
TP=cm2(1,1)
TN=cm2(2,2)
FP=cm2(1,2)
FN=cm2(2,1)

TPR=TP/(TP+FN)
FPR=FP/(FP+TN)
ACC=(TP+TN)/(TP+FN+FP+TN)
SPC=1-FPR
PPV=TP/(TP+FP)
NPV=TN/(FN+TN)



%% chi2 - kernel propio nuestro 
ker_train = 1-chi2_distance(train_featNorm,train_featNorm);
ker_test = 1-chi2_distance(test_featNorm,test_featNorm);
%train
model3 = svmtrain(yTrain(istrain,:), [(1:size(ker_train,1))' ker_train] , '-s 0 -t 4  -c 1 -w+1 1 ');
%predict
[predicted_label3, accuracy3, scores_train3] = svmpredict(yTrain(istrain,:), [(1:size(ker_train,1))' ker_train], model3);
[predicted_label3_t, accuracy3_t, scores_test3] = svmpredict(yTrain(istest,:), [(1:size(ker_test,1))' ker_test], model3);

%save tmp_data/chi2_scores predicted_label3 predicted_label3_t accuracy3 accuracy3_t scores_train3 scores_test3
tab2=ones(2,1);
tab = tabulate(yTrain(istest,1));
for i=1:size(tab,1)
    tab2(i,1)=tab(i,2);
end
cm=confusionmat(yTrain(istest,1),predicted_label3_t);
cm2=bsxfun(@rdivide,cm,tab2(:))*100

% Measures
TP=cm2(1,1)
TN=cm2(2,2)
FP=cm2(1,2)
FN=cm2(2,1)

TPR=TP/(TP+FN)
FPR=FP/(FP+TN)
ACC=(TP+TN)/(TP+FN+FP+TN)
SPC=1-FPR
PPV=TP/(TP+FP)
NPV=TN/(FN+TN)




%% RBF-chi2
gm = 1;
ker_train=exp(-gm*(1-ker_train));
ker_test=exp(-gm*(1-ker_test));
%train
model4 = svmtrain(yTrain(istrain,:), [(1:size(ker_train,1))' ker_train] , '-s 0 -t 4 -c 20 -w+1 1 ');
%predict
[predicted_label4, accuracy4, scores_train4]  = svmpredict(yTrain(istrain,:), [(1:size(ker_train,1))' ker_train], model4);
[predicted_label4_t, accuracy4_t, scores_test4] = svmpredict(yTrain(istest,:), [(1:size(ker_test,1))' ker_test], model4);

%save tmp_data/rbf-chi2_scores predicted_label4 predicted_label4_t accuracy4 accuracy4_t scores_train4 scores_test4

tab2=ones(2,1);
tab = tabulate(yTrain(istest,1));
for i=1:size(tab,1)
    tab2(i,1)=tab(i,2);
end
cm=confusionmat(yTrain(istest,1),predicted_label4_t);
cm2=bsxfun(@rdivide,cm,tab2(:))*100

% Measures
TP=cm2(1,1)
TN=cm2(2,2)
FP=cm2(1,2)
FN=cm2(2,1)

TPR=TP/(TP+FN)
FPR=FP/(FP+TN)
ACC=(TP+TN)/(TP+FN+FP+TN)
SPC=1-FPR
PPV=TP/(TP+FP)
NPV=TN/(FN+TN)

