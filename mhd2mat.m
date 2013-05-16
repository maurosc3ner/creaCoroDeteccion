%***************************
%*  mhd2mat
%*  
%***************************
%
%******************************

mkdir tmp_data;

clear all; clc; close all;
addpath src
addpath ../ReadData3D_version1k/mha
inDir = 'Training_seg/';

%% Flags
visualDebug=false;

MHD= dir(fullfile(inDir,'*.mhd'));
for j =1:numel(MHD),
    
    %% Files lecture

    refFilename=fullfile(inDir,[MHD(j).name(1:end-4) '.txt']);
    
    cprFilename=fullfile(inDir,MHD(j).name);
    reference=load(refFilename);
    info = mha_read_header(cprFilename)
    V = mha_read_volume(info);
    
    outputPath=fullfile(inDir,[MHD(j).name(1:end-4) '.mat']);
    
    save(outputPath,'reference','info','V')
    key=input('input key')
   
end
toc

