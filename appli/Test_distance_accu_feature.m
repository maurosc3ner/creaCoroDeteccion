%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Test distance accu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: E. Correa, june 02, 2013
% V: CPR.mhd
% Based on:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc; close all;
datasetDirectory='Training_SelVes/dt00'
addpath src

%% file lecture
refFilename=strcat(datasetDirectory,'/','vessel_5678.txt');
reference=load(refFilename);

% calling Ostium distance feature
dist=OstDistance(reference)
% attaching to the reference file
reference=[reference, dist]