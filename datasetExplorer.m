%*************************** 
%
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
inDir = 'Training_vessels/';

%busqueda de carpetas
DT= dir(fullfile(inDir,'dt*'));

%% Flags
visualDebug=false;
%% Training mode (decision criters)
% 0-> binary 0=healthy,soft|1=calc,mix
% 1-> binary 0=healthy,soft,mix|1=calc
% 2-> binary 0 healthy|1=calc,mix,narr>50
% 3-> multi 0=healthy|1=calc,mix|2=narr>50
% 4-> binary 0=healthy|1=grade narrowing>50
trainingMode=0;


yTrain=[];
numves=0;
svCounter=0;
spCounter=0;
spList=[];
svList=[];

L=5;
tic
%% dataset lecture
for j =1:numel(DT),
    
    %% Files lecture
    inDT=strcat(inDir,DT(j).name);
    %busqueda de vasos
    MHD= dir(fullfile(inDT,'*.mhd'));
    
    vesselTrainData=[];
    vesselyTrain=[];
    sickPatientFlag=0;
    for vessel_i=1:numel(MHD),
       
        refFilename=fullfile(inDT,[MHD(vessel_i).name(1:end-4) '.txt'])
        
        cprFilename=fullfile(inDT,MHD(vessel_i).name);
        reference=load(refFilename);
        info = mha_read_header(cprFilename)

        numves=numves+1;
        %Process
        [dims]=size(reference);
        %key=input('key');
        sliceOffset=round(dims(1)*.05);
        sickVesselFlag=0;
        for steps=sliceOffset:1:dims(1)-L-sliceOffset       % loop through CPR

            z0=round(L/2)+steps;
            % 3D Sampling pattern
            
            if trainingMode==0
                if reference(z0,6)>1.0
                    vesselyTrain=[vesselyTrain; 1];
                    sickVesselFlag=1;
                else
                    vesselyTrain=[vesselyTrain; 0];
                end
            elseif trainingMode==1
                if reference(z0,6)==2.0
                    vesselyTrain=[vesselyTrain; 1];
                else
                    vesselyTrain=[vesselyTrain; 0];
                end
                
            elseif trainingMode==2
                if reference(z0,6)>1.0
                    vesselyTrain=[vesselyTrain; 1];
                    sickVesselFlag=1;
                elseif reference(z0,7)>1
                    vesselyTrain=[vesselyTrain; 1];
                    sickVesselFlag=1;
                else
                    vesselyTrain=[vesselyTrain; 0];
                end
            elseif trainingMode==3
                if reference(z0,6)>1.0
                    vesselyTrain=[vesselyTrain; 1];
                elseif reference(z0,7)>1
                    vesselyTrain=[vesselyTrain; 2];
                else
                    vesselyTrain=[vesselyTrain; 0];
                end 
            elseif trainingMode==4
                if reference(z0,7)>1.0
                    vesselyTrain=[vesselyTrain; 1];
                    sickVesselFlag=1;
                else
                    vesselyTrain=[vesselyTrain; 0];
                end
            else
                vesselyTrain=[vesselyTrain; reference(z0,6)];
            end

            
        end   %rof process vessel
        
        if(sickVesselFlag==1)
            svCounter=svCounter+1;
            sickPatientFlag=1;
            svList{svCounter}={[DT(j).name MHD(vessel_i).name(1:end-4)]};
        end
        %key=input('key');
    end      %rof process Patient dataset
    if(sickPatientFlag==1)
            spCounter=spCounter+1;
            spList=[spList; DT(j).name];
    end
    
    %w = waitforbuttonpress;
    yTrain=[yTrain; vesselyTrain];
end
toc

svList
for i=1:numel(svList)
    svList{i}
end

[trainingMode numves svCounter spCounter]
