clear all;close all;clc;
addpath src
addpath ../ReadData3D_version1k/mha
inDir = 'Training_vessels/';

%busqueda de carpetas
DT= dir(fullfile(inDir,'dt*'));

tic;
for i =1:numel(DT),
    
    inDT=strcat(inDir,DT(i).name);
    %busqueda de vasos
    MHD= dir(fullfile(inDT,'*.mhd'));
    
    for j=1:numel(MHD),
       MHD(j).name 
       referenceFile=fullfile(inDT,[MHD(j).name(1:end-4) '.txt'])
       
    end
   
    
end
toc;