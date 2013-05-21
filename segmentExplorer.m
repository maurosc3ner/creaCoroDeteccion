%***************************
%*  3D Cylinder Model:
%*  I,Gr,Gt,Radon
%   Longitudinal 
%   Desgined for segment training
%*  2013,4,11
%*  References :
%*  Mittal10, 
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
%% Training mode (decision criters)
% 0-> binary 0=healthy,soft|1=calc,mix
% 1-> binary 0=healthy,soft,mix|1=calc
% 2-> binary 0 healthy|1=calc,mix,narr>50
% 3-> multi 0=healthy|1=calc,mix|2=narr>50
% 4-> binary 0=healthy|1=grade narrowing>50
trainingMode=4;

yTrain=[];
numseg=0;
ssCounter=0;
ssList=[];
L=5;
tic
%% Segment lecture
MHD= dir(fullfile(inDir,'*.mhd'));
for j =1:numel(MHD),
    
    %% Files lecture
    refFilename=fullfile(inDir,[MHD(j).name(1:end-4) '.txt'])
    reference=load(refFilename);
    [dims]=size(reference)
    numseg=numseg+1;
    %key=input('key input');
    
    sickSegFlag=0;
    segmentyTrain=[];
    for z0=round(L/2):1:dims(1)-round(L/2)+1     % loop through CPR

        if trainingMode==0
            if reference(z0,6)>1.0
                segmentyTrain=[segmentyTrain; 1];
                sickSegFlag=1;
            else
                segmentyTrain=[segmentyTrain; 0];
            end
        elseif trainingMode==1
            if reference(z0,6)==2.0
                segmentyTrain=[segmentyTrain; 1];
                sickSegFlag=1;
            else
                segmentyTrain=[segmentyTrain; 0];
            end
            
        elseif trainingMode==2
            if reference(z0,6)>1.0
                segmentyTrain=[segmentyTrain; 1];
                sickSegFlag=1;
            elseif reference(z0,7)>1
                segmentyTrain=[segmentyTrain; 1];
                sickSegFlag=1;
            else
                segmentyTrain=[segmentyTrain; 0];
            end
        elseif trainingMode==3
            if reference(z0,6)>1.0
                segmentyTrain=[segmentyTrain; 1];
                sickSegFlag=1;
            elseif reference(z0,7)>1
                segmentyTrain=[segmentyTrain; 2];
                sickSegFlag=1;
            else
                segmentyTrain=[segmentyTrain; 0];
            end
        elseif trainingMode==4
            if reference(z0,7)>1.0
                segmentyTrain=[segmentyTrain; 1];
                sickSegFlag=1;
            else
                segmentyTrain=[segmentyTrain; 0];
            end
        else
            segmentyTrain=[segmentyTrain; reference(z0,6)];
        end

    end   %rof process vessel
    
    if(sickSegFlag==1)
        ssCounter=ssCounter+1;
        ssList{ssCounter}={[MHD(j).name(1:end-4)]};
    end
    yTrain=[yTrain; segmentyTrain];
end
toc



hsList=[];
for i=1:numel(ssList)
    for j=1:numel(MHD)
        temp=char(ssList{i});
        if (strcmp(temp,MHD(j).name(1:end-4)))
            hsList=[hsList; j];
        end
    end
end

for (i=1:numel(hsList))
MHD(hsList(i)).name(1:end-4)
end

[trainingMode numseg ssCounter 14]
%     2   245    63    14