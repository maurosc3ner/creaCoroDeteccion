clear all; clc; close all;
mkdir models
%% Flags
%% Training mode 
% 0-> binary 0=healthy,soft|1=calc,mix
% 1-> binary 0=healthy,soft,mix|1=calc
trainingMode=3;

%% file lecture
if trainingMode==0
    trainFile='BC_L5_Train288'
    %trainFile='BC_Circle_allvesselsTrain'
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
else
    trainFile='MC_L5_allvesselsTrain'
end

load(strcat('tmp_data/',trainFile,'.mat'));

tabulate(yTrain)

%% Barajamiento de datos
part = cvpartition(yTrain,'holdout',0.4);
istrain = training(part); % data for fitting
istest = test(part); % data for quality assessment
tabulate(yTrain(istrain))

load 'models/rb500_TM3_seg_AP60'

cmpctRus=compact(rusTree)

sz(1)=whos('rusTree');
sz(2)=whos('cmpctRus');

[sz(1).bytes sz(2).bytes]

