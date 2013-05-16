%***************************
%*  filesaver for weka
%*  
%
%*  2013,4,11
%*  References :
%*  Chewla
%*  
%***************************
%
% Features (mean, std, gradMean, gradStd) are extracted with a 2D region of
% interest.
% ROI is defined by the circle equation:
%
%******************************


clear all; clc; close all;
addpath src

load 'tmp_data/BT_L5_allvesselsTrain.mat'



%dlmwrite('BT_L5_train.csv', [trainData logical(yTrain) ], 'delimiter', ',','newline','unix')
%M=csvread('BT_L5_train.csv',1,0);

%weka
matlab2csv('BT_L5_allvesselsTrain.csv',trainData,yTrain,2);

%addpath src/mweka
%arffwrite('BT_L5_train.arff','train',[],[],[trainData yTrain])


